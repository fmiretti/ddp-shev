function c = tensor_times_vector(a, b, dim)
%tensor_times_vector 
% 
% Parameters
% ----------
% a : double
%   A tensor.
% b : double
%   A vector.
% dim : numeric
%   dimension over which contraction is done.
% 
% Returns
% -------
% type
%   Explanation of anonymous return value of type ``type``.
% c : double
%   The product of a by b along the dimension dim.
%
% Examples
% --------
% C = tensor_times_vector(A, b, dim) 
%   multiplies the third-order tensor A by a column vector b along the 
%   dimension specified by dim.
%
% C = tensor_times_vector(A, b, dim) 
%   where A is a matrix and and b is a scalar performs matrix-times-vector
%   multiplication.
%
% C = tensor_times_vector(A, b, dim) 
%   where a and b are scalars performs scalar multiplication.

% Check that the vector size
if ~isequal(size(b), [size(a, dim) 1])
    error('The vector has a wrong size. Also, make sure it is a column vector.');
end

% Detect scalar case
if isscalar(a)
    c = a*b;
    return
end

% Get tensor size
sz = size(a);
% Fix for matrix-times-vector case
if ismatrix(a)
    sz = [sz 1];
end

% Permute the tensor so that the contracted dimension comes last
remdims = setdiff(1:numel(sz), dim);
a = permute(a, [remdims dim]);
sz = sz([remdims dim]);

% Reshape and multiply
a = reshape(a, prod(sz(1:end-1)), sz(end));
c = a * b;

% Convert the final result back to a matrix
c = reshape(c, sz(1:end-1));

end
