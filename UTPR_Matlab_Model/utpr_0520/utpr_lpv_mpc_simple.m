
clear
close all;

addpath(genpath('../common'))

%%
T  = 0.01;
nx = 6;
nu = 2;

x0 = [deg2rad(85); deg2rad(10); deg2rad(5); 0; 0; 0];
u_bar = [10; 10];   % 目前这个值确实太小了

Tf = 10; 

%%
N = 60;
Q = diag([1000, 2000, 1000, 10000, 200, 10]);
R = 300*diag([10, 9]);

%%
iter_max = 8; 
ref=[0*pi;0*pi;zeros(4,1)];

%% set MPC Parameters
% ref=[zeros(6,1)];
Q_hat = kron(eye(N), Q);
R_hat = kron(eye(N), R);

Acons = kron(eye(N), [eye(nu); -eye(nu)]);
bcons = kron(ones(2*N,1), u_bar);

% Initialize struct needed for dense problem formulation
Phi   = zeros(N*nx, nx);
Gamma = zeros(N*nx, N*nu);

%% Run simulation with MPC 

x_mpc = zeros(nx, 1);
u_mpc = zeros(nu, 1);
t_mpc = zeros(1, 1);
x_mpc(:,1) = x0;
x_hat(:,1) = ref - x0;

%% TODO: calculate com acceleration in x direction

rho = kron(ones(N+1,1), x0);
rho_conv = zeros(length(rho), iter_max+1);

cost_mpc = 0;

for j = 1:iter_max
    for i = 1:N
        % Calculate the system matrices of the LPV model at the current state and input
        x_range =  (i-1)*nx+1:i*nx;
        u_range =  (i-1)*nu+1:i*nu;
        [A, B] = utpr_lpv(rho(x_range), T);

        % Update state equation at the current time
        if i == 1
            Phi(x_range, :) = A;
        else
            Phi(x_range, :) = A*Phi(x_range - nx, :);
        end
        Gamma(x_range, u_range) = B;
        for l = 1:i-1
            z_range = (l-1)*nu+1:l*nu;
            Gamma(x_range, z_range) = A*Gamma(x_range-nx, z_range);
        end
    end

    % Calculate QP cost function matrices
    g = Gamma' * Q_hat * Phi * x_mpc(:,1);
    H = Gamma' * Q_hat * Gamma + R_hat;
    H = (H + H')/2;

    u = quadprog(H, g, Acons, bcons);

    if isempty(u)
        error('Quadprog failed to find a solution')
    end

    % Concatenate the current and predicted states
    rho = [x_mpc(:,1); Phi*x_mpc(:,1) + Gamma*u ]; % update rho
    rho_conv(:,j+1) = rho;
end



