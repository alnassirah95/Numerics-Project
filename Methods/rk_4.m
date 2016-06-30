function [ x_values, y_values ] = ...
    rk_4( f, initial_x, initial_y, step_size, steps )
% 4th order Runge-Kutta implementation for a single first-order
% differential equation.
%
% Syntax:
% [ xs, ys ] = rk_4( f, initial_x, initial_y, step_size, steps )
% - f: differential equation (some function handle)
% - initial_x: the value of x for the input initial condition.
% - initial_y: the value of y for the input initial condition.
% - step_size: the difference between x values for two consecutive steps
% - steps: the number of steps of RK-4 that the user wants to take
%
% Returns two matrices x_values and y_values where the nth row
% corresponds to the nth initial condition f(x_n)=y_n and stores the x and
% y values of the equation at each step.
%
% Example: Compute the differential equation f'(x) = x + y on the interval
% x = 0 to x = 1 using 10 steps each with a step size of 0.1. Use the
% initial condition f(0) = 0.
%
%   f_prime = @(x, y) x + y;
%   [xs, ys] = rk_4( f_prime, 0, 0, 0.1, 10 )

% By feeding in initial_x and initial_y as row vectors, we can compute
% RK-4 on multiple initial conditions at once.
% To do this, we first store the size of the initial_x and initial_y
% variables in the following two variables:
number_of_x_conditions = size(initial_x);
number_of_y_conditions = size(initial_y);

% Throughout RK-4 we will store x and y values in the following column
% vectors. We use steps + 1 because we also have to store the initial
% condition. The number of rows is equal to the number of initial
% conditions.
x_values = zeros(number_of_x_conditions(2), steps + 1);
y_values = zeros(number_of_y_conditions(2), steps + 1);

% The first value that goes in the column vectors is the initial condition.
x_values(1:end, 1) = initial_x;
y_values(1:end, 1) = initial_y;

% Store the current values of x and y
current_x = initial_x;
current_y = initial_y;

% Perform <steps> steps of RK-4 and store it in our column vectors.
% We go from 2 to steps + 1 so that we are always at the index where we
% will store the next x and y value in our column vectors x_vals, y_vals.
for step = 2:steps + 1
    % Calculate K1, K2, K3, and K4 using the formulas given by RK-4
    K1 = step_size * f(current_x, current_y);
    K2 = step_size * f(current_x + step_size / 2, current_y + K1 / 2);
    K3 = step_size * f(current_x + step_size / 2, current_y + K2 / 2);
    K4 = step_size * f(current_x + step_size, current_y + K3 );
    % Calculate new y value
    current_y = current_y + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
    % Calculate new x value
    current_x = current_x + step_size;
    % Store new x and y in column vectors
    x_values(1:end, step) = current_x;
    y_values(1:end, step) = current_y;
end

end

