function [ compressed_f ] = compress_function( f, start_x, end_x )
% Given a function f and an interval [start_x, end_x], returns a "squashed"
% version of f compressed_f such that all values of x in [start_x, end_x]
% have been compressed into the interval [-1, 1].
%
% Syntax:
% compressed_f = compress_function( f, start_x, end_x )
% - f: function handle for function we want to compress
% - start_x: value of x at the beginning of the interval
% - end_x: value of x at the end of the interval
%
% Note that this function works by creating a closure using the start_x and
% end_x variables, so the function that is created here should be manually
% cleared from memory.

% Create an error if start_x == end_x.
if start_x == end_x
    error( 'start_x == end_x; interval has length 0' )
end

% To compress, each value of x on [start_x, end_x] must map to x' via the
% equation x' = (2x - start_x - end_x) / (end_x - start_x), which implies
% that x = 1/2 ( x'(end_x - start_x) + start_x + end_x ). We use the latter
% equation to create the compressed function compressed_f.
compressed_f = @(x) f( 0.5 * (x * (end_x - start_x) + start_x + end_x) );

end

