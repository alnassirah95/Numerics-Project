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

% We store the current computed value of the integral in the following
% variable.
integral = 0;

% We store the equation f(z) / (z - a)^(n + 1) in an anonymous function so
% that we don't have to write the equation we are integrating over multiple
% times.
contour_fn = @(z) f(z) ./ ( (z - a) .^ (n + 1) );

% We store the angle around the point that we are currently at for in the
% following variable. This tells us where we should start integrating at
% the beginning of every batch later on.
current_angle = 0 + 0i;

% For extremely large numbers of steps, simply allocating memory for an
% offset for the beginning and end of each step can be extremely
% inefficient and may cause the program to crash. For that reason, we
% compute batchs of 10 million (1e7) steps at a time.
while steps > 0
    % Get the number of steps that should be computed in this batch. If the
    % number of remaining steps is greater than 10 million, we limit the
    % batch size to 10 million (1e7).
    batch_steps = steps;
    steps = steps - 1e7;
    if steps >= 1e7
        batch_steps = 1e7;
    end
    
    % Compute what the angle will be at the end of the current batch. If
    % this isn't our last batch (i.e. steps > 0) then we lower the angle.
    end_angle = 2 * pi;
    if steps > 0
        end_angle = 2 * pi / double( idivide( steps, uint16(10000000) ) );
    end
    
    % We store the offset from a for each step in the following variable. 
    % In order to compute the value of z at the beginning of each step, we
    % calculate a + offsets
    offsets = radius * exp( ...
        linspace( current_angle, end_angle, batch_steps + 1 ) * 1i );

    % We also store the offset from a for the end of each step in the 
    % following variable. Again, to calculate the value of z at the end of 
    % each step, we calculate a + end_offsets.
    end_offsets = offsets .* angle_shift;
    
    % Adjust the current angle for the next batch.
    current_angle = end_angle * angle_shift;
    
    % Now we compute the part of Simpson's rule that is "inside the
    % parentheses" for the current batch first and store it in a variable.
    batch_integral = contour_fn( offsets + a ) + ...
                    4 * contour_fn( half_angle_shift * offsets + a ) + ...
                    contour_fn( end_offsets + a );
                
    % Now we multiply by step size by subtracting offsets from end_offsets 
    % and taking the elementwise product.
    batch_integral = batch_integral .* (offsets - end_offsets);
    
    % Finally, we sum over all the steps. We put off dividing by 6 (per
    % Simpson's rule) until we have finished iterating for a little
    % additional efficiency.
    batch_integral = sum(batch_integral(:));
    
    % We update the integral variable with the value of the integral over
    % the current batch.
    integral = integral + batch_integral;
    
    % Since the offset and end_offset variables take up lots of memory for
    % high numbers of steps, we clear them before starting the next
    % iteration.
    clear offsets;
    clear end_offsets;
end
        
% Closures aren't automatically garbage-collected, so we have to clear the
% closure that we created (contour_fn) from memory.
clear contour_fn;

% We put off dividing by 6 (as according to Simpson's rule) while we were
% iterating; we do so now.
integral = integral / 6;

% Finally, now that we have computed the contour integral, we multiply by
% n! / (2 * pi * i) to get the nth derivative (by Cauchy's formula).
derivative = integral * factorial( n ) / (2 * pi * 1i);

end