function [ integral ] = trapezoidal_rule( f, start_x, end_x, steps )
% Uses the trapezoid rule to integrate a function over an interval.
%
% Syntax:
% integral = trapezoid_rule( f, start_x, end_x, steps )
% - f: function handle for function to integrate over
% - start_x: start of interval
% - end_x: end of interval
% - intervals: the number of steps to take to calculate the integral

% First, compute the step size for each iteration of the trapezoid rule
step_size = (end_x - start_x) / steps;

% Store the current computed value of the integral
integral = 0;

% Store the current value of x
current_x = start_x;

% Now we calculate the integral using the trapezoid rule.
for i = 1:steps
    % Calculate the new value of x at the end of the current step.
    step_x = current_x + step_size;
    
    % Update the integral variable using one step of the trapezoid rule
    integral = integral + step_size / 2 * ( f(current_x) + f(step_x) );
    
    % Update the value of x
    current_x = step_x;
end

end

