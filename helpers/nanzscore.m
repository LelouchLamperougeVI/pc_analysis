function Z = nanzscore(X, flag, dim)
% Stupid MATLAB can't zscore data with NaN's
% Same usage as zscore
%
% Function overloads:
%   Z = nanzscore(X)
%   Z = nanzscore(X, flag)
%   Z = nanzscore(X, 0, dim)
%   Z = nanzscore(X, flag, dim)
%
% Inputs:
%   X       Input matrix
%   flag    0: sample standard deviation (default)
%           1: population standard deviation
%   dim     dimension to operate on (default 1st non singleton)

switch nargin
    case 0
        error('No input was given');
    case 1
        flag = 0;
        dim = find( size(X) > 1, 1 );
    case 2
        dim = find( size(X) > 1, 1 );
    case 3
    otherwise
        error('too many input arguments');
end

Z = ( X - mean(X, dim, 'omitnan') ) ./ std( X, flag, dim, 'omitnan' );