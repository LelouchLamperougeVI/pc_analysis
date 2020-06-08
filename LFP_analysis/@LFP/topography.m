function topography(obj, FOV)
% does a bunch of stuff
% the important part is the obj.topo.loc which converts cell number indices
% to positional bin number

if nargin < 2
    FOV = obj.topo.FOV;
end

masks = arrayfun(@(x) obj.topo.maskNeurons == x, 1:size(obj.twop.deconv,2), 'uniformoutput', false);

obj.topo.centroid = cellfun(@regionprops, masks, 'uniformoutput',false);
obj.topo.centroid = arrayfun(@(x) reshape([obj.topo.centroid{x}.Centroid], 2, []), 1:length(obj.topo.centroid), 'uniformoutput',false);
obj.topo.centroid = cell2mat(obj.topo.centroid);

obj.topo.distances = distance(obj.topo.centroid);

[~,stack] = max( obj.analysis.stack );
obj.topo.loc = arrayfun(@(x) stack(x) .* (obj.topo.maskNeurons == x), 1:size(obj.twop.deconv,2), 'uniformoutput',false);
obj.topo.loc = sum( cell2mat( permute(obj.topo.loc, [1 3 2]) ), 3 );

    function d = distance(coor)
    %distances between centroids

    coor = coor .* FOV ./ size(obj.topo.mimg)';
    d = coor - permute(coor, [1 3 2]);
    d = sqrt( sum(d.^2, 1) );
    d = squeeze(d);

    end

end