function [A, B, M] = utpr_lpv(x, T)

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
% g = 9.81;


%% Disect states
% For easier calculation in the following lines, we split the states into theta and its derivative
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
%% Build preliminary matrices
% These are the matrices used to define also the nonlinear model in the paper. The vector G is not
% defined here, as it is replaced by the matrix K
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
% K matrix
% K = [-b1*si(theta(1))-b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b3*si(theta(1)+theta(2)+theta(3));
%      -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b2*si(theta(1)+theta(2))-b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b3*si(theta(1)+theta(2)+theta(3));
%      -b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b3*si(theta(1)+theta(2)+theta(3)), ...
%      -b3*si(theta(1)+theta(2)+theta(3)) ];

K = [-b1-b2-b3, ...
     -b2-b3, ...
     -b3;
     -b2-b3, ...
     -b2-b3, ...
     -b3;
     -b3, ...
     -b3, ...
     -b3 ];


% K = [-b1*si(theta(1))-b2*sin(theta(1)+theta(2))*(1/theta(1))-b3*sin(theta(1)+theta(2)+theta(3))*(1/theta(1)),   0,     0;...
%      0,   -b2*sin(theta(1)+theta(2))*(1/theta(2))-b3*sin(theta(1)+theta(2)+theta(3))*(1/theta(2)),      0;...
%      0,   0,      -b3*sin(theta(1)+theta(2)+theta(3))*(1/theta(3)) ];



% We cannot avoid inverting M, as it is needed for the input matrix B
Mi  = inv(M);
if any(~isfinite(Mi(:)))
    error('Hessian matrix Mi contains non-finite values.');
end

%% Assemble model matrices from parts
% From the previously defined parts, assemle the final model matrices A & B. At the same time, apply
% Euler's method for descritization.

A  = [  zeros(3)   eye(3) ;
       -Mi*K      -Mi*C   ];
B  = [ zeros(3,1)  zeros(3,1);
       Mi(:, 2)   Mi(:, 3)];
 
% Discretize if sampling time is given
if nargin >= 2
    A   = eye(size(A)) + T*A;
    B   = T*B;
end

if any(~isfinite(A(:)))
    error('Hessian matrix A contains non-finite values.');
end

if any(~isfinite(B(:)))
    error('Hessian matrix B contains non-finite values.');
end

if any(~isfinite(Mi(:)))
    error('Hessian matrix Mi contains non-finite values.');
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
