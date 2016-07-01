function [ row_vector ] = convert_to_row_vector( vector )
% Converts the input, if it is a column vector, into a row vector.
% Additionally, throws an exception if the input is not a row vector.
%
% Syntax:
% row_vector = convert_to_row_vector( initial_vector )

% Get the size of the vector
vector_dimensions = size( vector );

% Store the vector in row_vector
row_vector = vector;

% First check if the first dimension is greater than one
if vector_dimensions(1) ~= 1
    % Now check if the second dimension is greater than one. If it is, then
    % create an error.
    if vector_dimensions(2) ~= 1
        error( 'Bad input: not a vector' )
    % Otherwise, we transpose the row_vector variable so that now it is a
    % row vector.
    else
        row_vector = transpose(row_vector);
    end
end

end

