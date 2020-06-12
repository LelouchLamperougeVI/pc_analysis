function Y = bcircshift(A, k, dim, varargin)
% Stupic matlab `circshift' doesn't let you define a vector of shift values
% for each row/column
%
% Usage:
%   Y = bcircshift(A, k, dim, Name, Value)
%
% Inputs:
%   A:      input matrix (must be 2D)
%   k:      vector of shift values (must be same length as dimension being operated on)
%   dim:    dimension to operate on (default 1 - i.e. each column)
% Name, Value pairs:
%   'NaN'   'keep' (default) or 'ignore'
%           in case the matrix contains NaNs, choose whether to include
%           them with the shuffle or shuffle the values without touching
%           the NaNs
%
% Output:
%   Y:      shifted matrix

if ~ismatrix(A)
    error('input matrix must be 2D');
end

if nargin < 3
    dim = 1;
end

ops = parse_inputs(varargin);

if ops.nan
    nan_idx = any( isnan(A), mod(dim,2) + 1 );
    switch dim
        case 1
            A(nan_idx,:) = [];
        case 2
            A(:,nan_idx) = [];
    end
end

switch dim
    case 1
        Y = matcircshift(A, k);
    case 2
        Y = matcircshift(A', k)';
end

if ops.nan
    switch dim
        case 1
            tmp = NaN(length(nan_idx), size(Y, 2));
            tmp(~nan_idx, :) = Y;
        case 2
            tmp = NaN(size(Y,1), length(nan_idx));
            tmp(:, ~nan_idx) = Y;
    end
    Y = tmp;
end


function ops = parse_inputs(inputs)
ops.nan = false;

count = 1;
while(count < length(inputs))
    switch lower(inputs{count})
        case {'nan', 'nanflag'}
            switch lower(inputs{count+1})
                case {'none', 'keep'}
                    ops.nan = false;
                case {'ignore', 'omit', 'omitnan'}
                    ops.nan = true;
            end
    end
    count = count + 2;
end