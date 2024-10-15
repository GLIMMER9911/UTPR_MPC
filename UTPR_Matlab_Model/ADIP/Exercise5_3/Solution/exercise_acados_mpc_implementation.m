% Authors: Christian Hespe
% Institute of Control Systems, Hamburg University of Technology
% Date: 2021-11-16
%
% This file implements an MPC for a linearized model of the arm-driven
% inverted pendulum. It sets up the required quadratic program in
% single-shooting formulation manually.

addpath('..', genpath('../../common'))

clear
clear mex

setup_acados();

%% Model definition
dT = 0.01;
nx = 4; % Number of states
nu = 1; % Number of inputs

x0    = [deg2rad(45); 0; 0; 0];
u_bar = 1.6;

%% Controller Setup
ocp_N = 40;         % Length of the prediction horizon
ocp_T = ocp_N * dT; % Prediction horizon in seconds

iter_max = 10;

%% Define cost function
% The cost function on ACADOS is defined in terms of a virtual output y. In
% our case, we stack the states and inputs together to form y.

Q = diag([200, 1000, 0.1, 10]);
R = 2000;

ny = nx + nu;                   % Dimension of the virtual output
Vx = [ eye(nx); zeros(nu,nx) ]; % Mapping from states to virtual output
Vu = [ zeros(nx,nu); eye(nu) ]; % Mapping from inputs to virtual output

% Weighting matrix for the virtual output
W  = blkdiag(Q, R); 

%% Build ACADOS OCP model
% For details on the meaning of each of these terms, check the ACADOS documentation at 
% https://github.com/acados/acados/blob/master/docs/problem_formulation/problem_formulation_ocp_mex.pdf

ocp_model = acados_ocp_model();
ocp_model.set('name', 'adip');
ocp_model.set('T', ocp_T);

% Set dimension of the decision variables
ocp_model.set('dim_nx', nx);
ocp_model.set('dim_nu', nu);
ocp_model.set('dim_ny', ny);
ocp_model.set('dim_ny_e', 0); % No additional terminal constraint. See also the comment further down
ocp_model.set('dim_nz', 0);   % No algebraic variables

% Set symbolic variables
sym_x = casadi.SX.sym('x', nx);
sym_u = casadi.SX.sym('u', nu);
ocp_model.set('sym_x', sym_x);
ocp_model.set('sym_u', sym_u);

% Define dynamic properties of the model
ocp_model.set('dyn_type', 'explicit');    % Use explicit continuous-time ODE as model
ocp_model.set('dyn_expr_f', adip_ode(sym_x, sym_u));

% Set up cost function
ocp_model.set('cost_type', 'linear_ls');  % Linear least-squares cost function
ocp_model.set('cost_Vu', Vu);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vz', zeros(ny,0));    % No algebraic variables, so setting to 0
ocp_model.set('cost_W', W);
ocp_model.set('cost_y_ref', zeros(ny,1)); % Reference is at the origin, so setting to 0

% Set up input constraints
ocp_model.set('constr_C', zeros(1, nx));  % No constraints on the states
ocp_model.set('constr_D', 1);
ocp_model.set('constr_lg', -u_bar);       % Lower and upper limit on the control
ocp_model.set('constr_ug',  u_bar);

% Use no additional terminal constraint. Note that Acados already penalizes the state at the final
% step using the running cost. This is not immediately clear from the problem formulation.
ocp_model.set('cost_type_e', 'linear_ls');
ocp_model.set('cost_Vx_e', []);
ocp_model.set('cost_W_e', []);
ocp_model.set('cost_y_ref_e', []);

% Temporarily set initial conditions, without this the model contruction fails. Will be overwritten
% before the simulation anyway.
ocp_model.set('constr_x0', x0);

%% Configure ACADOS OCP solver
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', ocp_N);

% Set up the internal simulator. Here Euler's is selected, as this is the integrator used for qLMPC
% and CasADi, to make for a fair comparison.
ocp_opts.set('sim_method', 'erk');        % Use explicit Runge-Kutta ODE solver
ocp_opts.set('sim_method_num_stages', 4); % Use 4th order ODE solver

% Set up the NLP solver. If the solver is only allowed to perform a single iteration, switch to
% real-time iteration mode.
if iter_max == 1
    ocp_opts.set('nlp_solver', 'sqp_rti');
else
    ocp_opts.set('nlp_solver', 'sqp');
end
ocp_opts.set('nlp_solver_max_iter', iter_max); % Maximum number of iterations

% Set up the QP solver. Partial codensing is a mix of single- and
% multiple-shooting discretization of the problem in that multiplie start
% nodes are used, but they are each used for multiple prediction steps.
ocp_opts.set('qp_solver', 'partial_condensing_hpipm');
ocp_opts.set('qp_solver_cond_N', 10);

%% Generate OCP
ocp = acados_ocp(ocp_model, ocp_opts);

%% Run simulation with MPC
Tf   = 2;
Nsim = ceil(Tf / dT);

x_mpc = zeros(nx, Nsim+1);
u_mpc = zeros(nu, Nsim);
t_mpc = zeros(Nsim, 1);
x_mpc(:,1) = x0;

cost_mpc = 0;
for k = 1:Nsim
    tic
    ocp.set('constr_x0', x_mpc(:,k));
    ocp.solve();
    u_mpc(:,k) = clamp(ocp.get('u', 0), u_bar);
    t_mpc(k)   = toc;
    
    odefun = @(~, x) adip_ode(x,u_mpc(:,k));
    [~, x] = ode45(odefun, [(k-1)*dT, k*dT], x_mpc(:,k));
    x_mpc(:,k+1) = x(end,:);
    
    cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
end

fprintf('Total cost using MPC: %g\n', cost_mpc);
fprintf('Mean solver time :    %g [ms]\n', mean(t_mpc)*1e3);
fprintf('Median solver time :  %g [ms]\n', median(t_mpc)*1e3);
fprintf('Max solver time :     %g [ms]\n', max(t_mpc)*1e3);

%% Plot the results
figure(1)
stairs(0:dT:Tf-dT, u_mpc, 'r')
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
