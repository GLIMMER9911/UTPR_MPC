%% utpr CoM x test

clc;
clear;

syms m g b1 b2 b3 q1 q2 q3 qd1 qd2 qd3 qdd1 qdd2 qdd3;

cx = -1/(m*g)*(b1*sin(q1)+b2*sin(q1+q2)+b3*sin(q1+q2+q3));
cx_v = -1/(m*g)*(b1*cos(q1)*qd1 + b2*cos(q1+q2)*(qd1+qd2) + b3*cos(q1+q2+q3)*(qd1+qd2+qd3));

J_c = jacobian(cx_v, [qd1, qd2, qd3]);
cx_a = -1/(m*g)*(b1*cos(q1)*qdd1 + b2*cos(q1+q2)*(qdd1+qdd2) + b3*cos(q1+q2+q3)*(qdd1+qdd2+qdd3))+...
        1/(m*g)*(b1*sin(q1)*qd1*qd1 + b2*sin(q1+q2)*(qd1+qd2)*(qd1+qd2) + b3*sin(q1+q2+q3)*(qd1+qd2+qd3)*(qd1+qd2+qd3));

%Jcqd = 1/(m*g)*(b1*sin(q1)*qd1*qd1 + b2*sin(q1+q2)*(qd1+qd2)*(qd1+qd2) + b3*sin(q1+q2+q3)*(qd1+qd2+qd3)*(qd1+qd2+qd3))

Jcqd = jacobian((cx_a-J_c*[qdd1; qdd2; qdd3]), [qd1, qd2, qd3])