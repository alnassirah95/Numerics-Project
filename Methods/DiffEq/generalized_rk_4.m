function [ x_values, y_values ] = ...
    generalized_rk_4( function_array, initial_x, initial_ys, step_size, steps )
% Unlike the regular rk_4 method, this function is able to solve n
% differential equations simultaneously. However, this comes at the cost of
% ease of use, as this method must be fed a cell array of functions and
% each function must take in a vector as input.
%
% Syntax:
% xs, ys = generalized_rk_4( function_array, initial_xs, initial_ys, ...
%               step_size, steps )
% - function_array: a single-dimensional cell array of function handles
%       {u_1, u_2,..., u_n} of differential equations where the input to 
%       each equation is (x, y) where y is the vector given by
%
%           y = [u_1, u_2, ..., u_n]
%
%       That is, each function in the cell array is a function of x and a
%       vector containing the outputs of itself and all the other
%       functions.
% - initial_x: the value of x for our initial conditions.
% - initial_ys: the value of y for each equation at initial_x. The
%       first value corresponds to an initial condition for the first
%       function; the second value corresponds to an initial condition for
%       the second function; and so on.
% - step_size: the change in x between steps.
% - steps: the number of steps of RK-4 that we will take.
%
% Returns a vector of the x values at each step and an n by (<steps> + 1) 
% matrix, where the nth row corresponds to the solutions to the nth 
% differential equation.
%
% Example 1: RK-4 with three functions u1' = x + u2, u2' = x + u3, and 
% u3' = x + u1, with the initial conditions u1(0) = -1, u2(0) = 0, and 
% u3(0) = 1. Then, we plot the equations u1, u2, and u3.
% In this example, note that y is equivalent to the vector [u1, u2, u3].
% 
% u1_prime = @(x, y) x + y(2);
% u2_prime = @(x, y) x + y(3);
% u3_prime = @(x, y) x + y(1);
%
% initial_x = 0;
% initial_ys = [-1, 0, 1];
% step_size = 0.1;
% steps = 10;
%
% [xs, ys] = generalized_rk_4( {u1_prime, u2_prime, u3_prime}, ...
%                               initial_x, initial_ys, step_size, steps );
% plot(xs, ys)
%
% Example 2: here we generate 10 copies of the differential equation y' =
% x + y, each with different initial conditions, and plot them all on the
% same axes.
%
% initial_x = 0; initial_ys = linspace(0, 1, 10);
% steps = 10; step_size = 0.1;
%
% function_array = cell{10, 1};
% for i = 1:10
%   function_array{i} = @(x, y) x + y(i);
% end
% [xs, ys] = generalized_rk_4( function_array, initial_x, initial_ys, ...
%               step_size, steps );
% plot( xs, ys )

% Find the number of functions and the number of initial conditions.
num_functions = size( function_array );
num_conditions = size( initial_ys );

% Rearrange num_functions so that the number of functions is first, and the
% second value is (presumably) one (otherwise we will throw an error because
% we received a matrix of functions).
if num_functions(1) == 1
    num_functions = [num_functions(2), num_functions(1)];
end

% Rearrange num_conditions so that the number of initial conditions is
% first and the second value is 1 (if the user gave us a vector of initial
% conditions).
if num_conditions(1) == 1
    % We need the initial conditions to be a column vector, so we transpose
    % them (since the first dimension of the size matrix is 1).
    initial_ys = initial_ys';
    num_conditions = [num_conditions(2), num_conditions(1)];
end

% Error checking: check that the function_array and initial_conditions
% variables were given as vectors, and check that the number of functions
% and initial conditions are equal.
if num_functions(2) ~= 1
    error( 'Functions must be given as a vector.' )
elseif num_conditions(2) ~= 1
    disp( num_conditions )
    error( 'y values for initial conditions must be given as a vector.' )
elseif num_functions(1) ~= num_conditions(1)
    error( 'The number of initial conditions must equal the number of functions.' )
end

% Now that we have done all of our error-checking, we store the number of
% functions in num_functions (where we previously had the size matrix of
% the function_array variable).
num_functions = num_functions(1);

% Create a matrix that will store the value of x at each iteration.
x_values = linspace( initial_x, steps * step_size, steps + 1);

% We create a matrix to store the solutions to each of the equations. The
% nth row will store the solutions to the nth differential equation. We use
% steps + 1 rows because we also want to have a row for the initial
% conditions.
y_values = zeros( num_functions, steps + 1 );

% Set the first column of values in the y_values matrix equal to the
% initial conditions (which, from earlier, is a column vector).
y_values(:, 1) = initial_ys;

% Vectors to store the K values at each step.
K1 = zeros( num_functions, 1 );
K2 = zeros( num_functions, 1 );
K3 = zeros( num_functions, 1 );
K4 = zeros( num_functions, 1 );

% Now we iterate over the number of steps, computing the k values for each
% step.
for i = 1:steps
    % 1. Create K1 values.
    for j = 1:num_functions
        K1(j) = step_size * function_array{j}(x_values(i), ...
            y_values(:, i));
    end
    % 2. Create K2 values. First, store the vector of ys at which we
    % evaluate each function in the variable "k2_ys".
    k2_ys = y_values(:, i) + K1 / 2;
    for j = 1:num_functions
        K2(j) = step_size * function_array{j}(x_values(i) + step_size / 2, ...
            k2_ys );
    end
    % 3. Create K3 values. First, store the vector of ys at which we
    % evaluate each function in the variable "k3_ys".
    k3_ys = y_values(:, i) + K2 / 2;
    for j = 1:num_functions
        K3(j) = step_size * function_array{j}(x_values(i) + step_size / 2, ...
            k3_ys );
    end
    % 4. Create K4 values. First, store the vector of ys at which we
    % evaluate each function in the variable "k4_ys".
    k4_ys = y_values(:, i) + K3;
    for j = 1:num_functions
        K4(j) = step_size * function_array{j}(x_values(i) + step_size, ...
            k4_ys );
    end
    
    % 5. Use K1, K2, K3, and K4 values to compute the new values of y at
    % the current value of x.
    y_values(:, i + 1) = y_values(:, i) + (K1 + 2*K2 + 2*K3 + K4) / 6;
end

end

