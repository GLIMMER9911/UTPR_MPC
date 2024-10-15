%% UnderdrivernTriplePendulumMotor Parameters
PI = 3.141592653589793;
w = 2*PI;

% rand_para = rand()*0.1;

L1 = 0.325;
L2 = 0.2;
L3 = 0.2;
H = 0.05;
W = 0.05;

%% Link1 data
link1_m = 1.915;
link1_CM = [82.936*1e-3 0 0];
link1_MI = [0 0 25720.752*1e-6];
%link1_MI = [0 0 38896.101*1e-6];
% link1_MI = [8.333333333338782e-04 0.0271 0.0271];
link1_PI = [0 0 0];

%% Link2 data
link2_m = 1.469;
link2_CM = [65.468*1e-3 0 0];
link2_MI = [0 0 7667.511*1e-6];
%link2_MI = [0 0 13964.542*1e-6];
link2_PI = [0 0 0];

%% Link3 data
link3_m = 1.141;
link3_CM = [27.337*1e-3 0 0];
link3_MI = [0 0 4949.014*1e-6];
%link3_MI = [0 0 5801.831*1e-6];
link3_PI = [0 0 0 ];



















