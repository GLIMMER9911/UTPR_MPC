function u = utpr_mpc(x)
    %   UTPR_MPC
    %   UTPR 使用MPC进行控制
    %   x = (q1, q2, q3, dq1, dq2, dq3) state variable
    %   u = (torque1, torque2) input torque
    
    persistent ii N Q R iter_max ref Q_hat R_hat Acons bcons Phi Gamma rho t_mpc

    u = zeros(2,1);

    if isempty(ii)
        ii = 0;
    end   
    
    T  = 0.01;
    nx = 6;
    nu = 2;
    x0 = [deg2rad(0.6); deg2rad(0.6); deg2rad(0.6); 0; 0; 0];
    u_bar = [20; 20];   % 需要查文献确定取值

    %%
    N = 100;
    Q = diag([2000, 1000, 500, 2000, 100, 100]);
    R = 1*diag([10, 8]);

    %%
    iter_max = 2; 
    ref = [0*pi; 0*pi; zeros(4,1)];
    
    % ref=[zeros(6,1)];
    Q_hat = kron(eye(N), Q);
    R_hat = kron(eye(N), R);

    Acons = kron(eye(N), [eye(nu); -eye(nu)]);
    bcons = kron(ones(2*N,1), u_bar);
    
    % Initialize struct needed for dense problem formulation
    Phi   = zeros(N*nx, nx);
    Gamma = zeros(N*nx, N*nu);
    rho = kron(ones(N+1,1), x0);

    if ii == 0
        rho = kron(ones(N+1,1), x);
        %last_u = zeros(nu);
        t_mpc = 0;
    elseif ii > 0
        rho = kron(ones(N+1,1), x);
    end    
    
    tic
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
                Phi(x_range, :) = A*Phi(x_range-nx, :);
            end
            Gamma(x_range, u_range) = B;
            for l = 1:i-1
                z_range = (l-1)*nu+1:l*nu;
                Gamma(x_range, z_range) = A*Gamma(x_range-nx, z_range);
            end
        end
        
        max_value = 1e20;
        Gamma = max(Gamma, -max_value);  % 限制 Gamma 中的数值范围
        Gamma = min(Gamma,  max_value); 

        Phi = max(Phi, -max_value);  % 限制 Gamma 中的数值范围
        Phi = min(Phi,  max_value);

        regularization = 1e-10;
        Gamma = Gamma/(1 + regularization*norm(Gamma, "fro") );
        Phi = Phi/(1 + regularization*norm(Phi, "fro") );

        % Calculate QP cost function matrices
        g = Gamma' * Q_hat * Phi * x;
        H = Gamma' * Q_hat * Gamma + R_hat;
        H = (H + H')/2;
        
        options = optimoptions('quadprog', 'Algorithm', 'active-set');
        u = quadprog(H, g, Acons, bcons, [], [], [], [], [], [], options);
        
        if isempty(u)
            error('Quadprog failed to find a solution')
        end
        
        % Concatenate the current and predicted states
        rho = [ x; Phi*x + Gamma*u ]; % update rho 
    end
    toc
    t_mpc = toc + t_mpc;
    % u_mpc(:,k) = u(1:nu);

    % cost_mpc = cost_mpc + x_mpc(:,k)'*Q*x_mpc(:,k) + u_mpc(:,k)'*R*u_mpc(:,k);
    % rho = kron(ones(N+1,1), x_mpc(:,k+1));
    rho = [rho(nx+1:end); zeros(nx,1)];
    ii = ii + 1;
    u = clamp(u(1:nu), u_bar);

    % fprintf('Mean solver time :    %g [ms], step i: %g \n', (t_mpc/ii)*1e3, ii);

end

