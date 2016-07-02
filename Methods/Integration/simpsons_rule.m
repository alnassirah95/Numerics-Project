function [ integral ] = simpsons_rule( f, start_x, end_x, steps )
% Computes the integral over an interval using Simpson's rule.
%
% Syntax:
% integral = simpsons_rule( f, start_x, end_x, steps )
% - f: function handle for function to integrate over.
% - start_x: value of x at the beginning of the interval
% - end_x: value of x at the end of the interval
% - steps: number of steps of Simpson's rule to take within the interval

% Compute the step size for each step
step_size = (end_x - start_x) / steps;

% Store the current value of x in a variable
current_x = start_x;

% Store the current computed value of the integral in a variable
integral = 0;

% Iterate for the number of steps to compute the integral using Simpson's
% rule.
for i = 1:steps
    % Compute the value of x at the end of the current step
    end_x = current_x + step_size;
    
    % Compute the value of x in the middle of the current step
    middle_x = current_x + step_size / 2;
    
    % Now update the integral variable using Simpson's rule
    integral = integral + step_size / 6 * ...
        ( f(current_x) + 4 * f(middle_x) + f(end_x) );
    
    % Now update the current value of x to the value of x at the end of the
    % interval.
    current_x = end_x;
end

end

