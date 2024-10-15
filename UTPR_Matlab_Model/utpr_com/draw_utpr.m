function draw_utpr(state, spec)
%DRAW_ADIP Sketches an accurate representation of the ADIP for animation
%purposes.
%   Given the state of the ADIP at a single time instance, this function
%   generates an accurate visual representation.
%
%   Arguments:
%       state: Current state of the UTPR
%       spec:  (optional) Linespec (e.g. color) for the sketch

l_1  = 0.325;  % Length of the arm
l_2  =  0.2;   % Length of the pendulum
l_3  = 0.2;   % Length of the pendulum

%% Plotting the state
% Calculate the positions of the tips of arm and pendulum
x    = zeros(4,1);
y    = zeros(4,1);
x(2) = -l_1 * sin(state(1));
y(2) = l_1 * cos(state(1));
x(3) = x(2) - l_2 * sin(state(2));
y(3) = y(2) + l_2 * cos(state(2));
x(4) = x(3) - l_3 * sin(state(3));
y(4) = y(3) + l_3 * cos(state(3));

% Plot the lines representing arm and pendulum
if nargin <= 1
    plot(x, y, 'LineWidth', 1.5)
else
    plot(x, y, spec, 'LineWidth', 1.5)
end

% Add dots to represent the joints and boxes for the main housing
hold on
plot(x(1:3), y(1:3), 'k.', 'MarkerSize', 15, 'HandleVisibility', 'off')
hold off
rectangle('Position',[-0.05 -0.05 0.01 0.01])
rectangle('Position',[-0.25 -0.25 0.05 0.05])
end
