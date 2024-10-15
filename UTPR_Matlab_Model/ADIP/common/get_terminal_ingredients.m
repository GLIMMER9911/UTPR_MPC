function [W,P,F] = get_terminal_ingredients(A,B,Q,R,u_max,num_iter)
%GET_TERMINAL_CONSTRAINT_ITER Computes the three main components of an MPC 
%iteratevely, to guarantee stability. 
%A nominal F controller (LQ controller), A terminal
%state domain (LMI2), and a terminaFstate weight Lyapunov matrix.
%
%   This function requires a working copy of YALMIP in the Matlab path,
%   together with an appropriate solver for SDPs (e.g. SDPT-3 or Mosek).
%
%   This function will follow the same notation as in the thesis described
%   in [Loefberg, 2001]
%
%   Arguments:
%   N         Prediction Horizon
%   Q,R       Weighting Matrices
%   u_max     Maximum input control
%   Return values:
%   P         Quadratic terminaFstate weight
%   W         Condition for the terminal state domain constraint
%   F         Terminal controller

offset = 1e-8;

nx = size(A,1);
W  = [];
 

%% Step 1: Solve the Feasibility Problem (5.33) to get an initial feasible controller F:
% Feasibility Problem:
Y = sdpvar(nx,nx,'symmetric');
X = sdpvar(1,nx,'full');
t = sdpvar(1);


LMI_FP = [    Y     Y*A'+X'*B'       Y             X'     ;...
           A*Y+B*X      Y        zeros(nx)    zeros(nx,1) ;...
              Y      zeros(nx)     inv(Q)     zeros(nx,1) ;...
              X     zeros(1,nx)  zeros(1,nx)     inv(R)   ];

% Constraints settings: LMI_FP

Constraints = [ LMI_FP >= 0,... 
                   Y   >= t*eye(size(Y)),...
                   t   >=  offset];

% Solve optimization for LMI1_P
opts = sdpsettings('verbose', 0);
sol_LMI_P = optimize(Constraints, -t, opts);

if sol_LMI_P.problem ~= 0
    warning('Yalmip returned an error: "%s"', strip(yalmiperror(sol_LMI_P.problem)))
    F = [];
else
    
    F = value(X)/value(Y);
end

%% Step 2: Fix F and solve the terminal refion problem by maximizing det(Z):

Y = sdpvar(nx,nx,'symmetric');
Z = sdpvar(nx);
t = sdpvar(1);


% Terminal Region Problem

for index_iter = 1:num_iter

    LMI_Z1 = [     Z     Z*(A+B*F)' ;...
                (A+B*F)*Z     Z    ];

    LMI_Z2 = [  u_max^(2)-t  F*Z;...
                  Z*F'        Z];

    LMI_Z3 = [   Y       Y*(A+B*F)'      Y           Y*F'    ;...
              (A+B*F)*Y      Y        zeros(nx)   zeros(nx,1);...
                 Y        zeros(nx)     Q^(-1)    zeros(nx,1);...
                F*Y      zeros(1,nx)  zeros(1,nx)    R^(-1)];

    % STEP 2: Maximize the det(Z) or volume of the ellipsoid
    
    if rem(index_iter,2) == 1

        Constraints_X = [LMI_Z1 >= 0,...
                         LMI_Z2 >= 0,...
                         LMI_Z3 >= 0,...
                         t   >= offset,...
                         Y   >= offset*eye(size(Y)),...
                         Z   >= offset*eye(size(Z))];
        
        opts = sdpsettings('verbose', 0, 'solver', 'lmilab');
        sol_LMI_X = optimize(Constraints_X, -geomean(Z), opts);
        
        if sol_LMI_X.problem ~= 0
            
            fprintf('The last succesful iteration was %d.\n',index_iter-1)
            error('Yalmip returned an error: "%s"', strip(yalmiperror(sol_LMI_X.problem)))
        
        else
            
            if index_iter == num_iter
                
                Y = value(Y);
                Z = value(Z);
                break
                
            end
               
            F = sdpvar(1,nx,'full');    %Decision variable for next iter
            t = sdpvar(1);              %Value to maximinze for next iter
            Y = value(Y);               %Fixing resulting Y for next iter
            Z = value(Z);               %Fixing resulting Z for next iter

        end
           
    else
        
        %STEP 3: Maximize variable t and F
        
        Constraints_X = [LMI_Z1 >= 0,...
                         LMI_Z2 >= 0,...
                         LMI_Z3 >= 0,...
                            t   >=  offset];
        opts = sdpsettings('verbose', 0);
        sol_LMI_X = optimize(Constraints_X, -t, opts);
        
        if sol_LMI_X.problem ~= 0
            
            fprintf('The last succesful iteration was %d.\n',index_iter-1)
            error('Yalmip returned an error: "%s"', strip(yalmiperror(sol_LMI_X.problem)))
              
        else
            
            if index_iter == num_iter
            
                F = value(F);
                break
            
            end
            
            Y = sdpvar(nx,nx,'symmetric');    %Decision Variable next iter
            Z = sdpvar(nx);                   %Decision Variable next iter
            F = value(F);                     %Fix resulting F for next iter
            t = sdpvar(1);                    %Value to maximinze for next iter

        end
                


    end
    
end

  % Final Step: Obtain Original Variables:
  
  P = inv(Y);
  W = inv(Z);

end