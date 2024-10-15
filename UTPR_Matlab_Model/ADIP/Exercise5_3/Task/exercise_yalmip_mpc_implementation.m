% exercise_yalmip_mpc_implementation.m
%
% Authors: Christian Hespe & Santiago Pulido
% Institue of Control Systems, Hamburg University of Technology
% Email address: christian.hespe@tuhh.de
% Date: 2021-11-16
%
% This file presents the implementation of an MPC controller scheme,
% employing YALMIP toolbox from (Model Predictive Control Basics Example.)

% Solver ode45 for the differential equations of the pendulum was included.
% Note also that this file includes an animation of the behavior of the
% pendulum given the initial conditions and the parameters of the
% controller
%
% This script requires a working copy of YALMIP in the Matlab path,
% together with an appropriate solver for SDPs (e.g. SDPT-3 or Mosek).
%

addpath(genpath('../../common'))

clear

% Pendulum State Space Model Matrices:
T      = 0.01;
[A, B] = adip_lti(T);
nx     = size(A, 1);
nu     = size(B, 2);

x0    = [deg2rad(45); 0; 0; 0];
u_bar  = 1.6;

%% Explicit Prediction Form.
% MPC controller settings:
N = % FILL OUT
Q = % FILL OUT
R = % FILL OUT

% YALMIP variables definition:
U     = sdpvar(repmat(nu,1,N),ones(1,N));
x_hat = sdpvar(nx,1);

constraints = [];
objective   = 0;

% Objective and constraint functions generation:
x = x_hat;
for k = 1:N
    x = A*x + B*U{k};
    objective   = objective + % FILL OUT
    constraints = % FILL OUT
end

% Generate prepared solver for this MPC problem
opts = sdpsettings();
opts.quadprog.Algorithm = 'interior-point-convex';
controller = optimizer(constraints, objective, opts, x_hat, U{1});

%% Run simulation with MPC
Tf    = 2;
Nsim  = ceil(Tf / T);

x_mpc = zeros(nx, Nsim+1);
u_mpc = zeros(nu, Nsim);
t_mpc = zeros(Nsim,1);
x_mpc(:,1) = x0;

cost_mpc = 0;
for k = 1:Nsim
     tic
     [uk, err] = controller{x_mpc(:,k)};
     t_mpc(k)  = toc;
     
     if err ~= 0
        error('YALMIP can not find a solution')
     end
     u_mpc(:,k) = clamp(uk, u_bar);
     
     odefun = @(~, x) adip_ode(x,uk);
     [~, x] = ode45(odefun, [(k-1)*T, k*T], x_mpc(:,k));
     x_mpc(:,k+1) = x(end,:);
     
     cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
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
    draw_adip(x_mpc(:,k))
    
    grid on
    xlim([-0.1, 0.25])
    ylim([-0.1, 0.25])
    
    title('Movement of the Pendulum')
    xlabel('x Coordinate')
    ylabel('y Coordinate')
    
    drawnow limitrate
end
