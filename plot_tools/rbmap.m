function cmap = rbmap(varargin)
% colour the axis h in nice looking red-blue colormap
% 'caxis'       [min max] range - default clim from h
% 'normalize'   normalize the caxis at 0 - default true
% 'interp'      upsampling colourmap; must be odd - default 65 points
% 'colorbar'    add a colorbar - default true

[h, ops] = parse_inputs(varargin);

cmap = [    0.0196    0.1882    0.3804
            0.1294    0.4000    0.6745
            0.2627    0.5765    0.7647
            0.5725    0.7725    0.8706
            0.8196    0.8980    0.9412
            0.9686    0.9686    0.9686
            0.9922    0.8588    0.7804
            0.9569    0.6471    0.5098
            0.8392    0.3765    0.3020
            0.6980    0.0941    0.1686
            0.4039         0    0.1216    ];

cmap = interp1(1:size(cmap,1), cmap, linspace(1, size(cmap,1), ops.interp), 'spline');
cmap = (cmap - min(cmap(:))) ./ range(cmap(:));
idx = linspace(-1, 1, size(cmap,1));
if ops.normalize
    cmap = cmap( idx > ops.caxis(1) & idx < ops.caxis(2), : );
end

colormap(h, cmap)

if ops.colorbar
    p = get(h, 'Position');
    colorbar;
    set(h, 'Position', p);
end


function [h, ops] = parse_inputs(inputs)

ops.colorbar = true;
ops.caxis = [-1 1];
ops.normalize = true;
ops.interp = 65;

if isempty(inputs)
    h = gca;
    return
end

if isobject(inputs{1})
    h = inputs{1};
    idx = 2;
else
    h = gca;
    idx = 1;
end

while idx < length(inputs)
    switch lower(inputs{idx})
        case 'caxis'
            ops.caxis = inputs{idx+1};
        case 'normalize'
            ops.normalize = inputs{idx+1};
        case 'interp'
            ops.interp = inputs{idx+1};
        case 'colorbar'
            ops.colorbar = inputs{idx+1};
    end
    idx = idx + 2;
end

if ~mod(ops.interp, 2)
    warning('interp must be odd, adding 1 to given value');
    ops.interp = ops.interp + 1;
end

set(h, 'clim', ops.caxis);
