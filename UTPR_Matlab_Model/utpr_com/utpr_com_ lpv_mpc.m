% This file implements the qLMPC algorithm for the arm-driven inverted
% pendulum. The quadratic program is formulated in single-shooting form.

% TODO: add CoM acceleration constraint
% Jc*[qdd1, qdd2, qdd3] + Jcqd = desire_com_a;
% MPC求解的是关节的状态和力矩，如果想加的话，要在整个状态空间里面加


clear
close all;

addpath(genpath('../common'))

%%
T  = 0.01;
nx = 8;
nu = 2;

x0 = [deg2rad(85); deg2rad(10); deg2rad(5); 0; 0; 0; 0; 0];
u_bar = [10; 10];   % 目前这个值确实太小了


%%
N = 60;
Q = diag([1000, 2000, 1000, 10000, 200, 10, 0, 0]);
R = 300*diag([10, 9]);

%%
iter_max = 8; 
ref=[0*pi; 0*pi; zeros(6,1);];
% ref=[zeros(6,1)];
Q_hat = kron(eye(N), Q);
R_hat = kron(eye(N), R);

Acons = kron(eye(N), [eye(nu); -eye(nu)]);
bcons = kron(ones(2*N,1), u_bar);

% Initialize struct needed for dense problem formulation
Phi   = zeros(N*nx, nx);
Gamma = zeros(N*nx, N*nu);

%% Run simulation with MPC
Tf = 10;  
Nsim = ceil(Tf / T);

x_mpc = zeros(nx, Nsim+1);
u_mpc = zeros(nu, Nsim);
t_mpc = zeros(Nsim, 1);
x_mpc(:,1) = x0;
x_hat(:,1) = ref - x0;

%% TODO: calculate com acceleration in x direction

rho = kron(ones(N+1,1), x0);
rho_conv = zeros(length(rho), iter_max+1, Nsim);

cost_mpc = 0;
for k = 1:Nsim
    rho_conv(:,1,k) = rho;
    
    tic % 开启秒表计时
    for j = 1:iter_max
        for i = 1:N
            % Calculate the system matrices of the LPV model at the current state and input
            x_range =  (i-1)*nx+1:i*nx;
            u_range =  (i-1)*nu+1:i*nu;
            [A, B] = utpr_com(rho(x_range), T);

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

        u = quadprog(H, g, Acons, bcons);
        
        if isempty(u)
            error('Quadprog failed to find a solution')
        end
        
        % Concatenate the current and predicted states
        rho = [x_mpc(:,k); Phi*x_mpc(:,k) + Gamma*u ]; % update rho 
        rho_conv(:,j+1,k) = rho;
    end
    t_mpc(k)   = toc;
    u_mpc(:,k) = clamp(u(1:nu), u_bar);
    
    odefun = @(~, x) utpr_ode(x,u_mpc(:,k));
    [~, x] = ode45(odefun, [(k-1)*T, k*T], x_mpc(:,k));
    x_mpc(:,k+1) = unwrap(x(end,:));
    
    x_hat(:,k+1)=ref-x_mpc(:,k+1);
    % cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
    cost_mpc = cost_mpc + x_hat(:,k)'*Q*x_hat(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
    rho = [ rho(nx+1:end); zeros(nx,1) ];
end

fprintf('Total cost using MPC: %g \n', cost_mpc);
fprintf('Mean solver time :    %g [ms]\n', mean(t_mpc)*1e3);
fprintf('Median solver time :  %g [ms]\n', median(t_mpc)*1e3);
fprintf('Max solver time :     %g [ms]\n', max(t_mpc)*1e3);

%% Plot the results
figure(1)
stairs(0:T:Tf-T, u_mpc(1,:), 'r')
hold on
stairs(0:T:Tf-T, u_mpc(2,:), 'b')

yline( u_bar(1), 'k--')
yline(-u_bar(1), 'k--')
ylim(1.1*[-u_bar(1), u_bar(1)])

title('Control Trajectory')
xlabel('Time t')
ylabel('Torque \tau')
legend('u1','u2')
% figure(2)
% stairs(1:Nsim, t_mpc*1e3)
% title('MPC Solver Time')
% xlabel('Simulation step k')
% ylabel('Time t_{solv} in ms')

figure(2)
pause(0.1)
for k = 1:Nsim
    draw_utpr(x_mpc(:,k), 'r')
    
    xlim([-2, 2])
    ylim([-2, 2])
    
    title('Movement of the Pendulum')
    xlabel('x Coordinate')
    ylabel('y Coordinate')
    
    drawnow limitrate
end

figure(3)
subplot(131)
stairs(0:N-1, rho_conv(1:nx:N*nx,:,1))
xlabel('Timestep k')
ylabel('Arm angle \theta_1(k)')
subplot(132)
stairs(0:N-1, rho_conv(2:nx:N*nx,:,1))
xlabel('Timestep k')
ylabel('Pendulum angle \theta_2(k)')
subplot(133)
stairs(0:N-1, rho_conv(3:nx:N*nx,:,1))
xlabel('Timestep k')
ylabel('Pendulum angle \theta_3(k)')
legend(cellfun(@(c) ['Iteration ' num2str(c)], num2cell(0:iter_max), 'UniformOutput', false))
