% Authors: Christian Hespe
% Institute of Control Systems, Hamburg University of Technology
% Date: 2021-11-16
%
% This file implements an MPC for a linearized model of the arm-driven
% inverted pendulum. It sets up the required quadratic program in
% single-shooting formulation manually.

addpath(genpath('../../common'))

clear

%% Model definition
T  = 0.01;
[A, B] = adip_lti(T);
nx = size(A, 1);
nu = size(B, 2);

x0    = [deg2rad(45); 0; 0; 0];
u_bar = 1.6;

%% Controller Setup
N = 40;
Q = diag([200, 1000, 0.1, 10]);
R = 2000;

Q_bar = kron(eye(N), Q);
R_bar = kron(eye(N), R);

Phi   = zeros(N*nx, nx);
Gamma = zeros(N*nx, N*nu);
for k = 1:N
    x_range = (k-1)*nx+1:k*nx;

    Phi(x_range, :) = A^k;
    for j = 1:k
        Gamma(x_range, j) = A^(k-j)*B;
    end
end

H     = Gamma' * Q_bar * Gamma + R_bar;
g_hat = Gamma' * Q_bar * Phi;

Acons = kron(eye(N), [eye(nu); -eye(nu)]);
bcons = kron(ones(2*N,1), u_bar);

%% Run simulation with MPC
Tf   = 2;
Nsim = ceil(Tf / T);

x_mpc = zeros(nx, Nsim+1);
u_mpc = zeros(nu, Nsim);
t_mpc = zeros(Nsim, 1);
x_mpc(:,1) = x0;

cost_mpc = 0;
for k = 1:Nsim
    tic
    g = g_hat * x_mpc(:,k);
    u = quadprog(H, g, Acons, bcons);
    t_mpc(k) = toc;
    
    if isempty(u)
        error('Quadprog failed to find a solution')
    end
    u_mpc(:,k) = clamp(u(1:nu), u_bar);
    
    odefun = @(~, x) adip_ode(x,u_mpc(:,k));
    [~, x] = ode45(odefun, [(k-1)*T, k*T], x_mpc(:,k));
    x_mpc(:,k+1) = x(end,:);
    
    cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
end

fprintf('Total cost using MPC: %g\n', cost_mpc);
fprintf('Mean solver time :    %g [ms]\n', mean(t_mpc)*1e3);
fprintf('Median solver time :  %g [ms]\n', median(t_mpc)*1e3);
fprintf('Max solver time :     %g [ms]\n', max(t_mpc)*1e3);

%% Apply LQR to the linearized model
F = -dlqr(A, B, Q, R);

x_lqr = zeros(nx, Nsim+1);
u_lqr = zeros(nu, Nsim);
x_lqr(:,1) = x0;

cost_lqr = 0;
for k = 1:Nsim
    u_lqr(:,k) = clamp(F*x_lqr(:,k), u_bar);
    
    odefun = @(~, x) adip_ode(x,u_lqr(:,k));
    [~, x] = ode45(odefun, [(k-1)*T, k*T], x_lqr(:,k));
    x_lqr(:,k+1) = x(end,:);
    
    cost_lqr = cost_lqr + x_lqr(:,k)'*Q*x_lqr(:,k) + u_lqr(:,k)'*R*u_lqr(:,k);
end

fprintf('Total cost using LQR: %g\n', cost_lqr);

%% Plot the results
figure(1)
stairs(0:T:Tf-T, u_mpc, 'r')
hold on
stairs(0:T:Tf-T, u_lqr, 'b')
hold off

yline( u_bar, 'k--')
yline(-u_bar, 'k--')
ylim(1.1*[-u_bar, u_bar])

title('Control Trajectory')
xlabel('Time t')
ylabel('Torque \tau')
legend('MPC', 'LQR')

figure(2)
stairs(1:Nsim, t_mpc*1e3)
title('MPC Solver Time')
xlabel('Simulation step k')
ylabel('Time t_{solv} in ms')

figure(3)
pause(0.1)
for k = 1:Nsim+1
    draw_adip(x_mpc(:,k), 'r')
    hold on
    draw_adip(x_lqr(:,k), 'b')
    hold off
    
    grid on
    xlim([-0.1, 0.25])
    ylim([-0.1, 0.25])
    
    title('Movement of the Pendulum')
    xlabel('x Coordinate')
    ylabel('y Coordinate')
    legend('MPC', 'LQR')
    
    drawnow limitrate
end

figure(4)
spy(H)
title('Sparsity pattern of the Hessian')
