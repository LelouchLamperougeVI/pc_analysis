function h = violin(varargin)
% Make pretty violin plot
%
% h = violin(X)
% h = violin(X, 'nomedian', 'scatter', 'box', 'colours', c, 'labels', l)
%
% Input:
%   X:              cell array of data vectors
%
% Optional parameters:
%   'nomedian':     do not plot median and interquartile range lines
%   'scatter':      overlay scatter plot of data points
%   'box':          overlay boxplot
%   'colours', c:   custom colours specified as an Mx3 RGB-triplets matrix
%   'labels', l:    cell array of string labels
%
% by HaoRan Chang, PhD candidate
% Polaris Research Group, Canadian Centre for Behavioural Neuroscience
% University of Lethbridge, Alberta, Canada
[X, ops] = parse_ops(varargin);

if nargout > 0
    h = figure;
end
hold on

for ii = 1:length(X)
    [d, xi] = ksdensity(X{ii});
    d = (d - min(d)) ./ range(d) ./ 3;
    patch([d -d(end:-1:1)] + ii, [xi xi(end:-1:1)], ops.cmap(ii, :), 'facealpha', .4);
    
    if ops.line
        stats = [quantile(X{ii}, 1/4), median(X{ii}), quantile(X{ii}, 3/4)];
        idx = knnsearch(xi', stats');
        plot([d(idx(1)) -d(idx(1))] + ii, [xi(idx(1)) xi(idx(1))], ':', 'color', ops.cmap(ii, :), 'linewidth', 1.5);
        plot([d(idx(2)) -d(idx(2))] + ii, [xi(idx(2)) xi(idx(2))], '-', 'color', ops.cmap(ii, :), 'linewidth', 1.5);
        plot([d(idx(3)) -d(idx(3))] + ii, [xi(idx(3)) xi(idx(3))], ':', 'color', ops.cmap(ii, :), 'linewidth', 1.5);
    end
    
end
lims = get(gca, 'ylim');

linx = arrayfun(@(x) x .* ones(length(X{x}), 1), 1:length(X), 'uniformoutput', false);
try
    linx = cell2mat(linx);
catch
    linx = cell2mat(linx');
end
try
    liny = cell2mat(X);
catch
    liny = cell2mat(X');
end

if ops.scatter
    scatter(linx(:), liny(:), 32, repelem(ops.cmap, cellfun(@length, X), 1), 'filled', 'jitter', 'on', 'jitteramount', .15, 'markerfacealpha', .5);
end
if ops.box
    boxplot(liny(:), linx(:), 'plotstyle', 'compact', 'colors', ops.cmap);
end

xticks(1:length(X));
xticklabels(ops.labels);
xlim([0 length(X)+1]);
ylim(lims);


function [X,ops] = parse_ops(inputs)
X = inputs{1};
if ~isa(X, 'cell') && isa(X, 'double')
    X = {X};
end
ops.line = true;
ops.scatter = false;
ops.box = false;
ops.cmap = distinguishable_colors(length(X), {'w', 'k'});
ops.labels = strsplit(num2str(1:length(X)));

count = 2;
while count <= length(inputs)
    switch lower(inputs{count})
        case {'nomed', 'nomedian'}
            ops.line = false;
        case 'scatter'
            ops.scatter = true;
        case {'box', 'boxplot'}
            ops.box = true;
        case {'colors', 'colours'}
            ops.cmap = inputs{count + 1};
            count = count + 1;
        case 'labels'
            ops.labels = inputs{count + 1};
            count = count + 1;
        otherwise
            error(['''' inputs{count} ''' is not a valid option.']);
    end
    count = count + 1;
end

if ops.box
    ops.line = false;
end

if ops.box && ops.scatter
    warning('Boxplots and scatter should not be used together.');
end