function p = padj(p, type)
% Adjust p-values for multiple comparisons by Bonferroni or other methods.
%
% Inputs
%   p: vector of p-values. If matrix, operates on the first dimension.
%   type: adjustment method (default - 'bonferroni')
%
% Output
%   p: adjusted p-values
%
% List of methods:
%   'bonferroni'
%   'holm-bonferroni'
%   'hochberg'
%   'sidak'
%   'tch' ()
%   'benjamini–hochberg'
%
%
% Implementation of methods described in:
%
% Blakesley, R. E., Mazumdar, S., Dew, M. A., Houck, P. R., Tang, G., Reynolds, C. F., 3rd, & Butters, M. A. (2009).
% Comparisons of methods for multiple hypothesis testing in neuropsychological research.
% Neuropsychology, 23(2), 255–264.
% https://doi.org/10.1037/a0012850
%
% Copyright (C) 2022 - HaoRan Chang
% SPDX-License-Identifier: GPL-2.0-or-later


if nargin < 2
    type = 'bonferroni';
end

switch lower(type)
    case {'bon', 'bonferroni'}
        m = size(p, 1);
        p = m .* p;
        
    case {'holm', 'holm-bon', 'holm-bonferroni'}
        m = (size(p, 1) : -1 : 1)';
        p = step_through(p, m, 'down');
        
    case {'hoch', 'hochberg'}
        m = (1 : size(p, 1))';
        p = step_through(p, m, 'up');
        
    case {'sid', 'sidak'}
        m = size(p, 1);
        p = 1 - (1 - p) .^ m;
        
    case {'tch', 'tukey-ciminera-heyse'}
        m = size(p, 1);
        p = 1 - (1 - p) .^ sqrt(m);
        
    case {'bh', 'b-h', 'fdr', 'ben', 'benjamini', 'benjamini–hochberg'}
        m = size(p, 1) ./ (size(p, 1) : -1 : 1)';
        p = step_through(p, m, 'up');
        
    case {'none'}
    otherwise
        error('Invalid adjustment method.');
end

p = cat(3, p, ones(size(p)));
p = min(p, [], 3);


function p = step_through(p, m, direction)
% step-up, step-down, step-sideways, you name it, we'll do it :)

switch lower(direction)
    case 'down'
        [pp, order] = sort(p, 'ascend');
        fun = @max;
    case 'up'
        [pp, order] = sort(p, 'descend');
        fun = @min;
end

pp = m .* pp;

for ii = 1:size(pp, 1)
    pp(ii, :) = fun(pp(1:ii, :), [], 1);
end

[~, order] = sort(order);
cols = repmat(1:size(pp, 2), size(pp, 1), 1);
idx = sub2ind(size(pp), order, cols);
p = pp(idx);
