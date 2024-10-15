function [A, B] = adip2_lpv(x, T)
%ADIP_LTI Linearized ADIP model
%   This function provides system matrices for the ADIP that are obtained
%   by Jacobi linearization around the origin of the state-space.
%
%   Arguments:
%       T: (optional) Sampling time for discrete-time model
%           A continuous-time model is returned if T is not given

%% Define constants
l_1  = 0.325;  % Length of the arm
lg_1 =  0.245436; % Center of gravity of the arm
l_2  =  0.4;   % Length of the pendulum
lg_2 = 0.292805;
m_1  = 1.915;    % Mass of the arm
m_h  = 130e-3;    % Mass of the encoder
m_2  =  1.61;  % Mass of the pendulum 2.6
g    = 9.81;      % Force of gravity

fv_1 = 1.3e-3;    % Coefficient of friction for the arm
fv_2 = 2.2e-5;    % Coefficient of frcition for the pendulum

J1 = 0.025720752;
J2 = 32871.952e-6;

theta = x(1:2);
omega = x(3:4);

%% Define intermediate matrices
m21 = m_2*l_1*l_2;
M   = [ (m_1+m_2)*l_1^2             m21*cos(theta(1)-theta(2)) ;
        m21*cos(theta(1)-theta(2))  m_2*l_2^2                  ];
C   = [  fv_1                                 m21*sin(theta(1)-theta(2))*omega(2) ;
        -m21*sin(theta(1)-theta(2))*omega(1)  fv_2                                ];
K   = [ -(m_1*lg_1+(m_h+m_2)*l_1)*g*si(theta(1))   0                        ;
         0                                          -m_2*l_2*g*si(theta(2)) ];

Mi  = inv(M);

%% Define state-space model
A  = [  zeros(2)   eye(2) ;
       -Mi*K      -Mi*C   ];
B  = [ zeros(2,1) ; Mi(:, 2) ];

% Discretize if sampling time is given
if nargin >= 1
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