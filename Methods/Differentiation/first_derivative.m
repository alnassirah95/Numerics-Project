function [ derivative ] = first_derivative( f, x, h )
% Computes the first derivative of an input function f using the five point
% midpoint formula:
%
%   f'(x) = 1 / 12h ( f(x-2h) - 8f(x - h) + 8f(x + h) - f(x + 2h) )
%
% The error using this method is O(h^4).
%
% Syntax:
% derivative = first_derivative( f, x, h )
% - f: function handle for function to find derivative of.
% - x: point at which derivative is to be evaluated.
% - h: "Step size" used to compute the derivative; smaller values of h are
%      more accurate (but can become less numerically stable).

% Use formula listed in description in order to compute the first
% derivative:
derivative = (f(x - 2*h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / (12 * h);

end

