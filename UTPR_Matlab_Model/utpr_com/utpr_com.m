function [A, B] = utpr_com(x, T)
%%  CoM dynamics calculate
%   x(q1, q2, q3, \dot{q1}, \dot{q2}, \dot{q3}, cx, cx_v ) 
%   dx(\dot{q1}, \dot{q2}, \dot{q3}, \ddot{q1}, \ddot{q2}, \ddot{q3}, cx_v, cx_a)

g = 9.81; % Force of gravity
l_1 = 0.325;  % Length of the arm
lg_1 = 0.245436; % Center of gravity of the arm
l_2 = 0.2; % Length of the pendulum
lg_2 = 0.165468;
l_3 = 0.2; % Length of the pendulum
lg_3 = 0.127337;
m1 = 1.915; % Mass of the arm
m2 = 1.469; % Mass of the pendulum
m3 = 1.141; % Mass of the pendulum
J1 = 0.025720752;
J2 = 0.007667511;
J3 = 0.004949014;

%% 
theta = x(1:3);
omega = x(4:6);

%% intermediate variable
a1 = J1+m1*lg_1^2+(m2+m3)*l_1^2;
a2 = J2+m2*lg_2^2+m3*l_2^2;
a3 = (m2*lg_2+m3*l_2)*l_1;
a4 = J3+m3*lg_3^2;
a5 = m3*l_1*lg_3;
a6 = m3*l_2*lg_3;

b1 = (m1*lg_1+m2*l_1+m3*l_1)*g;
b2 = (m2*lg_2+m3*l_2)*g;
b3 = (m3*lg_3)*g;

%% Equations of motion
M  = [ a1+a2+a4+2*a5*cos(theta(2)+theta(3))+2*a3*cos(theta(2))+2*a6*cos(theta(3)), ... 
        a2+a4+a3*cos(theta(2))+a5*cos(theta(2)+theta(3))+2*a6*cos(theta(3)), ...
        a4+a5*cos(theta(2)+theta(3))+a6*cos(theta(3));
        a2+a4+a3*cos(theta(2))+a5*cos(theta(2)+theta(3))+2*a6*cos(theta(3)), ...
        a2+a4+2*a6*cos(theta(3)), ...
        a4+a6*cos(theta(3));
        a4+a5*cos(theta(2)+theta(3))+a6*cos(theta(3)), ...
        a4+a6*cos(theta(3)), ...
        a4  ];

C  = [  -a5*(omega(2)+omega(3))*sin(theta(2)+theta(3))-a3*omega(2)*sin(theta(2))-a6*omega(3)*sin(theta(3)), ...
        -a5*(omega(1)+omega(2)+omega(3))*sin(theta(2)+theta(3))-a3*(omega(1)+omega(2))*sin(theta(2))-a6*omega(3)*sin(theta(3)), ...
        -(omega(1)+omega(2)+omega(3))*(a5*sin(theta(2)+theta(3))+a6*sin(theta(3)));
        a5*omega(1)*sin(theta(2)+theta(3))+a3*omega(1)*sin(theta(2))-a6*omega(3)*sin(theta(3)), ...
        -a6*omega(3)*sin(theta(3)), ...
        -a6*(omega(1)+omega(2)+omega(3))*sin(theta(3));
        a5*omega(1)*sin(theta(2)+theta(3))+a6*(omega(1)+omega(2))*sin(theta(3)), ...
        a6*(omega(1)+omega(2))*sin(theta(3)), ...
        0  ];

G   = [ -b1*sin(theta(1))-b2*sin(theta(1)+theta(2))-b3*sin(theta(1)+theta(2)+theta(3));
        -b2*sin(theta(1)+theta(2))-b3*sin(theta(1)+theta(2)+theta(3));
        -b3*sin(theta(1)+theta(2)+theta(3))];

K = [-b1*si(theta(1))-b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
     -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
     -b3*si(theta(1)+theta(2)+theta(3));
     -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
     -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
     -b3*si(theta(1)+theta(2)+theta(3));
     -b3*si(theta(1)+theta(2)+theta(3)), ...
     -b3*si(theta(1)+theta(2)+theta(3)), ...
     -b3*si(theta(1)+theta(2)+theta(3)) ];

%% CoM dynamics
% cx = -1/(m*g)*(b1*sin(theta(1))+b2*sin(theta(1)+theta(2))+b3*sin(theta(1)+theta(2)+theta(3)));

J_c = [-(b2*cos(theta(1) + theta(2)) + b1*cos(theta(1)) + b3*cos(theta(1) + theta(2) + theta(3)))/(g*m), ...
        -(b2*cos(theta(1) + theta(2)) + b3*cos(theta(1) + theta(2) + theta(3)))/(g*m), ...
        -(b3*cos(theta(1) + theta(2) + theta(3)))/(g*m)];

J_cv = [(2*b1*omega(1)*sin(theta(1)) + 2*b2*sin(theta(1) + theta(2))*(omega(1) + omega(2)) + 2*b3*sin(theta(1) + theta(2) + theta(3))*(omega(1) + omega(2) + omega(3)))/(g*m), ... 
    (2*b2*sin(theta(1) + theta(2))*(omega(1) + omega(2)) + 2*b3*sin(theta(1) + theta(2) + theta(3))*(omega(1) + omega(2) + omega(3)))/(g*m), ...
    (2*b3*sin(theta(1) + theta(2) + theta(3))*(omega(1) + omega(2) + omega(3)))/(g*m)];

% cx_a = J_c / M *[0;u] + (- J_c / M * C + J_cv)*omega - J_c / M * K * theta;

%% Assemble model matrices from parts
Mi = inv(M);

A  = [  zeros(3)       eye(3)           zeros(3,2) ;
       -Mi*K          -Mi*C             zeros(3,2) ;
        zeros(1,3)     J_c              zeros(1,2) ;
        -J_c/M*K     -J_c/M*C+J_cv      zeros(1,2) ];
B  = [ zeros(3,1)      zeros(3,1);
       Mi(:, 2)        Mi(:, 3);
       zeros(1, 1)     zeros(1, 1)
       J_c*Mi(:,2)     J_c*Mi(:,3)];

% Discretize if sampling time is given
if nargin >= 2
    A   = eye(size(A)) + T*A;
    B   = T*B;
end
end

function y = si(x)
%SI Defines the SI function
%   Matlab does not define the SI function, so we need to. Just taking sin(x)/x will result in a
%   singularity for x = 0, which is undesirable

    if x == 0
        y = 1;
    else
        y = sin(x) / x;
    end
end


end

