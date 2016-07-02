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

% The angle by which we should rotate around the point at which we are
% trying to evaluate the derivative is equal to 2 * pi / steps. Since we
% are integrating counterclockwise, we take the conjugate of this value.
angle_shift = conj( cos(2 * pi / steps) + sin(2 * pi / steps) * 1i );

% Half of the above angle, used to calculate the offset at the middle of
% each step. Again, we take the conjugate of this angle.
half_angle_shift = conj( cos(pi / steps) + sin(pi / steps) * 1i );

% We store the offset from a for each step in the following variable. In
% order to compute the value of z at the beginning of each step, we
% calculate a + offsets
offsets = radius * exp( linspace( 0, 2 * pi, steps + 1 ) * 1i );

% We also store the offset from a for the end of each step in the following
% variable. Again, to calculate the value of z at the end of each step, we
% calculate a + end_offsets.
end_offsets = offsets .* angle_shift;

% We store the equation f(z) / (z - a)^(n + 1) in an anonymous function so
% that we don't have to write the above equation multiple times.
contour_fn = @(z) f(z) ./ ( (z - a) .^ (n + 1) );

% Now we compute the part of Simpson's rule that is "inside the
% parentheses" first and store it in a variable.
integral = contour_fn( offsets + a ) + ...
            4 * contour_fn( half_angle_shift * offsets + a ) + ...
            contour_fn( end_offsets + a );

% Now we multiply by step size by subtracting offsets from end_offsets and
% taking the elementwise product.
integral = integral .* (offsets - end_offsets);

% Finally, we sum over all the steps and divide by 6 (due to Simpson's 
% rule) in order to finish computing the integral over the circle.
integral = sum(integral(:)) / 6;
        
% Closures aren't automatically garbage-collected, so we have to clear the
% closure that we created (contour_fn) from memory. We also clear offsets
% and end_offsets from memory since they take up a large amount of space
% for large numbers of steps.
clear contour_fn;
clear offsets;
clear end_offsets;

% Finally, now that we have computed the contour integral, we multiply by
% n! / (2 * pi * i) to get the nth derivative (by Cauchy's formula).
derivative = integral * factorial( n ) / (2 * pi * 1i);

end