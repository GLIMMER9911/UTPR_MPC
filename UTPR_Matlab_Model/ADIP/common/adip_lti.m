function [A, B, C] = adip_lti(T)
%ADIP_LTI Linearized ADIP model
%   This function provides system matrices for the ADIP that are obtained
%   by Jacobi linearization around the origin of the state-space.
%
%   Arguments:
%       T: (optional) Sampling time for discrete-time model
%           A continuous-time model is returned if T is not given

%% Define constants
l_1  = 13.95e-2;  % Length of the arm
lg_1 =  6.975e-2; % Center of gravity of the arm
l_2  =  7.8e-2;   % Length of the pendulum
m_1  = 115e-3;    % Mass of the arm
m_h  = 130e-3;    % Mass of the encoder
m_2  =  73.1e-3;  % Mass of the pendulum
fv_1 = 1.3e-3;    % Coefficient of friction for the arm
fv_2 = 2.2e-5;    % Coefficient of frcition for the pendulum
g    = 9.81;      % Force of gravity

%% Define intermediate matrices
J = [  l_2^2    -l_1*l_2          ;
      -l_1*l_2  (m_1/m_2+1)*l_1^2 ] / (m_1*l_1^2*l_2^2);
K = [ -(m_1*lg_1+(m_h+m_2)*l_1)*g   0         ;
       0                           -m_2*l_2*g ];
L = diag([fv_1, fv_2]);

%% Define state-space model
A = [  zeros(2)   eye(2) ;
       -J*K       -J*L    ];
B = [ zeros(2,1) ;
       J * [1; 0] ];
C = [ eye(2)  zeros(2) ];

% Discretize if sampling time is given
if nargin >= 1
    A   = eye(size(A)) + T*A;
    B   = T*B;
end
end
