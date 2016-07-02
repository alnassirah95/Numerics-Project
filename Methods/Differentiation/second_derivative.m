function [ derivative ] = second_derivative( f, x, h )
% Computes the second derivative using the three point midpoint formula,
% which is as follows:
%
%  f''(x) = 1 / h^2 (f(x - h) - 2f(x) + f(x + h))
%
% The error on this formula is O(h^2).
% Syntax:
% derivative = second_derivative( f, x, h )
% - f: function handle for function to find derivative of.
% - x: point at which derivative is to be evaluated.
% - h: "Step size" used to compute the derivative; smaller values of h are
%      more accurate (but can become less numerically stable).

% Compute the second derivative using the three point midpoint formula from
% the description.
derivative = (f(x - h) - 2*f(x) + f(x + h)) / (h.^2);

end

