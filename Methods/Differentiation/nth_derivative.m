function [ derivative ] = nth_derivative( f, a, n, steps, radius )
% Computes the nth derivative using the Cauchy integral formula. This
% method should only be used if no other formula exists for computing an
% nth derivative. Returns a complex value.
%
% We integrate the given function f / (z-x)^(n+1) over a circle of the 
% input radius in the complex plane using Simpson's rule and multiply by 
% n! / (2 * pi * i); this gives us the derivative at the input point.
% 
%
% Syntax:
% derivative = nth_derivative( f, a, n, h )
% - f: function handle for function to find derivative of.
% - a: point at which derivative is to be evaluated.
% - n: the number of the derivative that we should compute.
% - steps: the number of steps that should be used to evaluate the
%      integral.
% - radius: radius of the circle in the complex plane that we will
%      integrate over.
% Both the radius and steps variables affect the accuracy of the integral;
% too many steps or too small a radius may cause numerical errors.

% Error checking: 
% 1. Check that the number of the derivative is a positive integer.
% 2. Check that the steps variable is a positive integer.
% 3. Check that the radius is positive.
if n <= 0
    error( 'Number n of derivative must be a positive integer.' )
elseif steps <= 0
    error( 'Number of steps must be a positive integer.' )
elseif radius <= 0
    error( 'Radius must be positive.' )
end

% In order to compute the integral, we store a current offset from the
% point at which we are trying to compute the integral, which starts off
% being equal to radius. Then at each iteration we shift z by multiplying
% it by another variable, angle_shift.
offset = radius;

% The angle by which we should rotate around the point at which we are
% trying to evaluate the derivative is equal to 2 * pi / steps. Since we
% are integrating counterclockwise, we take the conjugate of this value.
angle_shift = conj( cos(2 * pi / steps) + sin(2 * pi / steps) * 1i );

% Half of the above angle, used to calculate the offset at the middle of
% each step. Again, we take the conjugate of this angle.
half_angle_shift = conj( cos(pi / steps) + sin(pi / steps) * 1i );

% Store the value of the integral on the contour
integral = 0;

% We store the equation f(z) / (z - a)^(n + 1) in an anonymous function so
% that we don't have to write the above equation over and over again.
contour_fn = @(z) f(z) / ( (z - a).^(n + 1) );

% Now we iterate over the number of steps in order to compute the integral.
for i = 1:steps
    % Get the offset from the point at which we are trying to differentiate
    % that is at the end of the current step.
    end_offset = offset * angle_shift;
    % Get the offset from the point at which we are trying to differentiate
    % that is at the middle of the current step.
    middle_offset = offset * half_angle_shift;
    % Use Simpson's rule to approximate the integral between offset and
    % end_offset. Our "step size" here is equal to (a + end_offset) - (a +
    % offset) = end_offset - offset, so we multiply by that. We put off
    % dividing by 6 until the end for increased efficiency.
    integral = integral + (end_offset - offset) * ...
        (contour_fn(a + offset) + 4 * contour_fn(a + middle_offset) + ...
         contour_fn(a + end_offset));
    % Update the current offset before the next step.
    offset = end_offset;
end

% Closures aren't automatically garbage-collected, so we have to clear the
% closure that we created (contour_fn) from memory.
clear contour_fn;

% While we were integrating, we didn't divide by 6 in order to increase
% efficiency; we do so now.
integral = integral / 6;

% Finally, now that we have computed the contour integral, we multiply by
% n! / (2 * pi * i) to get the nth derivative (by Cauchy's formula).
derivative = integral * factorial( n ) / (2 * pi * 1i);

end

