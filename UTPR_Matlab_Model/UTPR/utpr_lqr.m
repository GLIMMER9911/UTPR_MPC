% Authors: Donghao Li
% Institute of Control Systems, SUSTech
% Date: 2024-5-7
%
% This file implements the qLMPC algorithm for the arm-driven inverted
% pendulum. The quadratic program is formulated in single-shooting form.

clc;
clear;

addpath(genpath('./common'))

%%
T  = 0.01;
x0 = [deg2rad(0.5); deg2rad(0.5); deg2rad(0.5); 0; 0; 0];
[A, B] = utpr_lti(T, x0);
nx = 6;
nu = 2;

u_bar = [100; 100];   % 需要查文献确定取值，写成两个或三个

% %%
N = 50;
Q = diag([1000, 500, 200, 10000, 200, 100]);
R = diag([10, 10]);
iter_max = 4; 

Q_hat = kron(eye(N), Q);
R_hat = kron(eye(N), R);

Tf = 10;  
Nsim = ceil(Tf / T);

F = -dlqr(A, B, Q, R);
x_lqr = zeros(nx, Nsim+1);
u_lqr = zeros(nu, Nsim);
x_lqr(:,1) = x0;

ee_pos = zeros(2, Nsim + 1 );
ee_pos(:, 1) = utpr_xy_position(x0);

ref=[0.; 0.0; -0.0; zeros(3,1)];
cost_lqr = 0;
for k = 1:Nsim
    u_lqr(:,k) = clamp(F*(x_lqr(:,k)-ref), u_bar);

    odefun = @(~, x) utpr_ode(x,u_lqr(:,k));
    [~, x] = ode45(odefun, [(k-1)*T, k*T], x_lqr(:,k));
    x_lqr(:,k+1) = x(end,:);
    
    % a = -0.003;
    % b = 0.003;
    % A = [(a + ( b - a ) * rand(1,1));0;0; zeros(3,1);];
    % 
    % if (400 <= k) && ( k <= 480)
    %     x_lqr(:,k+1) = x_lqr(:,k+1) + A;
    % end
    cost_lqr = cost_lqr + x_lqr(:,k)'*Q*x_lqr(:,k) + u_lqr(:,k)'*R*u_lqr(:,k);

    ee_pos(:, k+1) = utpr_xy_position(x_lqr(:,k+1));

end

fprintf('Total cost using LQR: %g\n', cost_lqr);

figure(1)
for k = 1:Nsim
    draw_utpr(x_lqr(:,k), 'r')

    xlim([-1, 1])
    ylim([-0.5, 1])

    title('Movement of the Pendulum')
    xlabel('x Coordinate')
    ylabel('y Coordinate')
    pause(0.01)

end

 fillcolor1=[189, 30, 30]/256;
 fillcolor2=[252, 170, 103]/256;
 fillcolor3=[0, 70, 222]/256;

figure(2)
stairs( 0:T:Tf-T, u_lqr(1,:), Color=fillcolor1, LineWidth=1.6)
hold on
stairs( 0:T:Tf-T, u_lqr(2,:),  Color=fillcolor3, LineWidth=1.6)
hold off

% yline( u_bar(1), 'k--')
% yline(-u_bar(2), 'k--')
ylim([-4, 4])
title('Control Trajectory')
xlabel('Time t')
ylabel('Torque \tau')
legend('torque 2', 'torque 3')

figure(3)
stairs( 0:T:Tf, x_lqr(1,:), Color=fillcolor1, LineWidth=1.6)
hold on
stairs( 0:T:Tf, x_lqr(2,:), Color=fillcolor2, LineWidth=1.6)
hold on
stairs( 0:T:Tf, x_lqr(3,:), Color=fillcolor3, LineWidth=1.6)
hold off

title('Control Trajectory')
xlabel('Time t')
ylabel('Joint position')
legend('joint 1', 'joint 2', 'joint 3')


figure(4)
stairs( 0:T:Tf, x_lqr(4,:), Color=fillcolor1, LineWidth=1.6)
hold on
stairs( 0:T:Tf, x_lqr(5,:), Color=fillcolor2, LineWidth=1.6)
hold on
stairs( 0:T:Tf, x_lqr(6,:), Color=fillcolor3, LineWidth=1.6)
hold off

title('Control Trajectory')
xlabel('Time t')
ylabel('Joint Velocity')
legend('joint 1', 'joint 2', 'joint 3')

%% 将代码保存为.xlsx
% ee_pos = ee_pos';
% datapos = [x_lqr(1:3,:)', ee_pos];
% xlsdata = ["joint_pos1", "joint_pos2", "joint_pos3", "x", "y"; datapos ];
% xlswrite("E:\Code\MATLAB\UTPR_Matlab_Model\UTPR\LQR_result\UTPR_lqr_result.xlsx", xlsdata);




