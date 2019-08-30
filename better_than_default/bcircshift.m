function Y = bcircshift(A, k, dim)
% Stupic matlab `circshift' doesn't let you define a vector of shift values
% for each row/column
%
% Usage:
%   Y = bcircshift(A, k, dim)
%
% Inputs:
%   A:      input matrix (must be 2D)
%   k:      vector of shift values (must be same length as dimension being operated on)
%   dim:    dimension to operate on (default 1 - i.e. each column)
%
% Output:
%   Y:      shifted matrix

if ~ismatrix(A)
    error('input matrix must be 2D');
end

if nargin < 3
    dim = 1;
end

switch dim
    case 1
        Y = matcircshift(A, k);
    case 2
        Y = matcircshift(A', k)';
end