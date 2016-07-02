function [ integral ] = gaussian_quadrature( f, start_x, end_x )
% Use third-degree Gaussian quadrature to compute an integral.
%
% Syntax:
% integral = gaussian_quadrature( f, start_x, end_x )
% - f: function handle for function to integrate over.
% - start_x: value of x at the beginning of the interval.
% - end_x: value of x at the end of the interval.

% First, we have to "compress" f so that the interval [start_x, end_x]
% becomes the interval [-1, 1] for a new function compressed_f.
compressed_f = compress_function( f, start_x, end_x );

% Now we use the Gaussian quadrature to compute the integral.
% The third-degree Gaussian quadrature gives us the following way to
% compute the integral over the given interval:
% 
% integral = 8/9 f(-sqrt(3/5)) + 5/9 f(0) + 8/9 f(sqrt(3/5))

integral = ( end_x - start_x ) / 2 * ( 5/9 * compressed_f( -sqrt( 3/5 ) ) + ...
    8/9 * compressed_f(0) + 5/9 * compressed_f( sqrt( 3/5 ) ) );

% Since compress_function creates a closure using the variables start_x and
% end_x, we have to clear it from memory (as it is not automatically
% garbage collected).
clear compressed_f;

end

