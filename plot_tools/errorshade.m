function h = errorshade(varargin)
% plots nicer errorbar plots
%
% overloads:
%   errorshade(y, err, name-pairs)
%   errorshade(x, y, err, name-pairs)
%
% name,value pairs:
%   color/colour    rgb triplets and letters
%   line/linespec   e.g. '--r'
%   h/target        h
%   alpha           .5 default

[x, y, err, ops] = parse_inputs(varargin);

if isempty(ops.h)
    figure;
    h = gca;
else
    h = ops.h;
end
axes(h);
hold on
for ii = 1:size(y, 2)
    draw(h, x, y(:,ii), err(:,ii), ops.colour(ii,:), ops.alpha, ops.linespec);
end


function draw(h, x, y, err, colour, alpha, linespec)
fill(h, [x; x(end:-1:1)], [y+err; y(end:-1:1)-err(end:-1:1)], colour, 'EdgeColor',colour, 'FaceAlpha', alpha);
plot(h, x, y, 'color', colour, 'linestyle',linespec);


function [x, y, err, ops] = parse_inputs(inputs)

ops.colour = 'k';
ops.linespec = '-';
ops.h = [];
ops.alpha = .5;

track = 4;

idx = cellfun(@ischar, inputs);
idx = find(idx, 1);
if isempty(idx)
    idx = length(inputs);
    track = 3;
end

if idx < track
    y = inputs{1};
    err = inputs{2};
    x = 1:size(y,1);
else
    x = inputs{1};
    y = inputs{2};
    err = inputs{3};
end

if ~ismatrix(x) || ~ismatrix(y) || ~ismatrix(y)
    error('inputs not matrix')
end
if size(ops.colour,2)>1; ops.colour=ops.colour'; end
if size(x,2)>1; x=x'; end

if size(y,2) > 1
    ops.colour = distinguishable_colors(size(y,2));
end

while idx < length(inputs)
    switch lower(inputs{idx})
        case {'color', 'colour'}
            ops.colour = inputs{idx+1};
        case {'line', 'linespec'}
            ops.linespec = inputs{idx+1};
        case {'h', 'target'}
            ops.h = inputs{idx+1};
        case 'alpha'
            ops.alpha = inputs{idx+1};
        otherwise
            error([inputs{idx} ' is not a valid option']);
    end
    idx = idx + 2;
end