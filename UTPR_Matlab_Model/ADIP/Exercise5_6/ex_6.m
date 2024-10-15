% Authors: Christian Hespe
% Institute of Control Systems, Hamburg University of Technology
% Date: 2021-11-16
%
% This file implements the qLMPC algorithm for the arm-driven inverted
% pendulum. The quadratic program is formulated in single-shooting form.

addpath(genpath('../common'))

clear

%%
T  = 0.01;
nx = 4;
nu = 1;

x0    = [deg2rad(45); deg2rad(30); 0; 0];
u_bar = 1.6;

%%
N = 40;
Q = diag([200, 1000, 0.1, 10]);
R = 2000;
iter_max = 5;          %

Q_hat = kron(eye(N), Q);
R_hat = kron(eye(N), R);

Acons = kron(eye(N), [eye(nu); -eye(nu)]);
bcons = kron(ones(2*N,1), u_bar);

% Initialize struct needed for dense problem formulation
Phi   = zeros(N*nx, nx);
Gamma = zeros(N*nx, N*nu);

%% Run simulation with MPC
Tf   = 2;
Nsim = ceil(Tf / T);

x_mpc = zeros(nx, Nsim+1);
u_mpc = zeros(nu, Nsim);
t_mpc = zeros(Nsim, 1);
x_mpc(:,1) = x0;

rho = kron(ones(N+1,1), x0);
rho_conv = zeros(length(rho), iter_max+1, Nsim);

cost_mpc = 0;
for k = 1:Nsim
    rho_conv(:,1,k) = rho;
    
    tic
    for j = 1:iter_max
        for i = 1:N
            % Calculate the system matrices of the LPV model at the current state and input
            x_range =  (i-1)*nx+1:i*nx;
            u_range =  (i-1)*nu+1:i*nu;
            [A, B] = adip_lpv(rho(x_range), T);

            % Update state equation at the current time
            if i == 1
                Phi(x_range, :) = A;
            else
                Phi(x_range, :) = A*Phi(x_range-nx, :);
            end
            Gamma(x_range, u_range) = B;
            for l = 1:i-1
                z_range = (l-1)*nu+1:l*nu;
                Gamma(x_range, z_range) = A*Gamma(x_range-nx, z_range);
            end
        end

        % Calculate QP cost function matrices
        g = Gamma' * Q_hat * Phi * x_mpc(:,k);
        H = Gamma' * Q_hat * Gamma + R_hat;
        H = (H + H')/2;

        u   = quadprog(H, g, Acons, bcons);
        if isempty(u)
            error('Quadprog failed to find a solution')
        end
        
        % Concatenate the current and predicted states
        rho = [ x_mpc(:,k); Phi*x_mpc(:,k) + Gamma*u ]; % update rho 
        rho_conv(:,j+1,k) = rho;
    end
    t_mpc(k)   = toc;
    u_mpc(:,k) = clamp(u(1:nu), u_bar);
    
    odefun = @(~, x) adip_ode(x,u_mpc(:,k));
    [~, x] = ode45(odefun, [(k-1)*T, k*T], x_mpc(:,k));
    x_mpc(:,k+1) = x(end,:);
    
    cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
    rho = [ rho(nx+1:end); zeros(nx,1) ];
end

fprintf('Total cost using MPC: %g\n', cost_mpc);
fprintf('Mean solver time :    %g [ms]\n', mean(t_mpc)*1e3);
fprintf('Median solver time :  %g [ms]\n', median(t_mpc)*1e3);
fprintf('Max solver time :     %g [ms]\n', max(t_mpc)*1e3);

%% Plot the results
figure(1)
stairs(0:T:Tf-T, u_mpc, 'r')

yline( u_bar, 'k--')
yline(-u_bar, 'k--')
ylim(1.1*[-u_bar, u_bar])

title('Control Trajectory')
xlabel('Time t')
ylabel('Torque \tau')

figure(2)
stairs(1:Nsim, t_mpc*1e3)
title('MPC Solver Time')
xlabel('Simulation step k')
ylabel('Time t_{solv} in ms')

figure(3)
pause(0.1)
for k = 1:Nsim
    draw_adip(x_mpc(:,k), 'r')
    
    xlim([-0.25, 0.25])
    ylim([-0.25, 0.25])
    
    title('Movement of the Pendulum')
    xlabel('x Coordinate')
    ylabel('y Coordinate')
    
    drawnow limitrate
end

figure(4)
subplot(121)
stairs(0:N-1, rho_conv(1:nx:N*nx,:,1))
xlabel('Timestep k')
ylabel('Arm angle \theta_1(k)')
subplot(122)
stairs(0:N-1, rho_conv(2:nx:N*nx,:,1))
xlabel('Timestep k')
ylabel('Pendulum angle \theta_2(k)')
legend(cellfun(@(c) ['Iteration ' num2str(c)], num2cell(0:iter_max), 'UniformOutput', false))
