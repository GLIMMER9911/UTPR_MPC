function [W,P,F] = get_terminal_ingredients_lqr(A,B,Q,R,u_bar)
%GET_TERMINAL_INGREDIENTS_LQR Computes terminal ingredients using an LQ
%controller as the terminal controller.
%   This function requires a working copy of YALMIP in the Matlab path.
%
%   This function will follow the same notation as in the thesis described
%   in [Loefberg, 2001]
%
%   Arguments:
%   A, B      Model Matrices
%   Q, R      Weighting Matrices
%   u_max     Maximum control input
%   Return values:
%   P         Quadratic terminaFstate weight
%   W         Condition for the terminaFstate domain constraint

% Define offset that makes the LMIs strict
offset = 1e-8;

nx = size(A,1);
W  = [];

%% Step 1: Find the nominal LQ controller:
% First, use Matlab build-ins to solve the DARE. The resulting controller
% will work as a terminal controller and the solution P to the Riccati
% equation satisfies the LMI (5.28) from the lecture notes by equality.

[F, P] = dlqr(A,B,Q,R);

if isempty(F)
    warning('There is no solution to the LQR problem for this plant')
    return
else
    F = -F; % We use a different sign convention than Matlab does
end

%% Step 2: Terminal State Domain X LMI_X:
% Given the terminal penalty and terminal controller, find the ellipsoid
% that maximises the volume of the terminal domain.

Z   = sdpvar(nx);
LMI = [     Z      Z*(A+B*F)' ;...
        (A+B*F)*Z      Z      ];
Constraints = [ Z      >= offset*eye(nx) ,...
                LMI    >= 0              ,...
                F*Z*F' <= u_bar^2        ];

% Solve optimization for LMI_X
opts = sdpsettings('verbose', 0, 'solver', 'lmilab');
sol  = optimize(Constraints, -geomean(Z), opts);

if sol.problem ~= 0
    warning('Yalmip returned an error: "%s"', strip(yalmiperror(sol.problem)))
else
    W = inv(value(Z));
end
end
