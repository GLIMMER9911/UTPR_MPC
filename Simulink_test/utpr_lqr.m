function u = utpr_lqr(x)
    %   utpr_mpc
    %   utpr 使用 mpc 进行控制
    %   x = (q1, q2, q3, dq1, dq2, dq3) state variable
    %   u = (torque1, torque2) input torque

    persistent ii Q R
    u = zeros(2,1);

    if isempty(ii)
        ii = 0;
    end 

    T  = 0.01;
    % x0 = [deg2rad(0.4); deg2rad(0.4); deg2rad(0.4); 0; 0; 0];
    
    [A, B] = utpr_lti(T,x);
    nx = 6;
    nu = 2;

    u_bar = [20; 20];   

    N = 50;
    Q = diag([1000, 500, 200, 10000, 200, 100]);
    R = diag([10, 10]);

    F = -dlqr(A, B, Q, R);
    
    ref=[0.; 0.0; -0.0; zeros(3,1)];

    ii = ii+1;

    u = clamp(F*(x - ref), u_bar);

end

