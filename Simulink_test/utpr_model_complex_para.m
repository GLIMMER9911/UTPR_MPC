%% UnderdrivernTriplePendulumMotor Parameters
PI = 3.141592653589793;
w = 2*PI;

L1 = 0.325;
L2 = 0.2;
L3 = 0.2;
H = 0.05;
W = 0.05;

rand_para = (2*rand(6,1)-ones(6,1))*0.15;

%% Link1 data
link1_m = 1.915 * (rand_para(1,1) + 1);
link1_CM = [245.436 * 1e-3 0 0];
link1_MI = [0 0 25720.752*1e-6* (rand_para(2,1) + 1)];
link1_PI = [0 0 0];

%% Link2 data
link2_m = 1.469 * (rand_para(3,1) + 1);
link2_CM = [165.468 * 1e-3 0 0];
link2_MI = [0 0 7667.511*1e-6* (rand_para(4,1) + 1)];
link2_PI = [0 0 0];

%% Link3 data
link3_m = 1.141 * (rand_para(5,1) + 1);
link3_CM = [127.337*1e-3 0 0];
link3_MI = [0 0 4949.014*1e-6* (rand_para(6,1) + 1)];
link3_PI = [0 0 0 ];

%% 形状定义
% 参数定义
n = 50;
theta = linspace(pi/2, -pi/2, n);  % 角度从 0 到 π（半圆）

R_outer = 0.02;
R_inner = 0.01;

% 计算外圆和内圆的坐标
x_outer = R_outer * cos(theta); % 外圆 X 坐标
y_outer = R_outer * sin(theta); % 外圆 Y 坐标

x_inner = R_inner * cos(theta); % 内圆 X 坐标
y_inner = R_inner * sin(theta); % 内圆 Y 坐标

% 闭合路径的连接，顺时针连接外圆点，然后逆时针连接内圆点
left_end_x_section = [flipud(x_outer') flipud(y_outer'); x_inner' y_inner'];

% 第一根杆件
theta = linspace(-pi/2, pi/2, n);
points_x1 = R_inner*cos(theta) ;
points_y1 = R_inner*sin(theta) ;
points_x2 = R_outer*cos(theta) + (L1);
points_y2 = R_outer*sin(theta) ;
% link1_section = [0,0; -L1/2, 0; points_x1' points_y1'; -L1/2, 0.2; L1/2, 0.2; flipud(points_x2') flipud(points_y2')];
link1_section = [ 0, -R_outer; points_x2' points_y2'; 0 R_outer; flipud(points_x1') flipud(points_y1');];

% 第二根杆件
theta = linspace(-pi/2, pi/2, n);
points_x1 = R_inner*cos(theta) ;
points_y1 = R_inner*sin(theta) ;
points_x2 = R_outer*cos(theta) + (L2);
points_y2 = R_outer*sin(theta) ;
% link1_section = [0,0; -L1/2, 0; points_x1' points_y1'; -L1/2, 0.2; L1/2, 0.2; flipud(points_x2') flipud(points_y2')];
link2_section = [ 0, -R_outer; points_x2' points_y2'; 0 R_outer; flipud(points_x1') flipud(points_y1');];

link3_section = link2_section;


% % 检查形状
% figure;
% plot(link2_section(:,1), link2_section(:,2), '-or');
% axis equal;
% title('Cross-Section Shape');











