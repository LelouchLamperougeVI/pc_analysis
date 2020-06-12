function topography(obj)
% does a bunch of stuff
% the important part is the obj.topo.loc which converts cell number indices
% to positional bin number

mask = arrayfun(@(x) regionprops(any(obj.topo.maskNeurons == x, 3)), 1:size(obj.twop.deconv, 2), 'uniformoutput', false);
area = cellfun(@(x) max([x.Area]), mask);
centroids = arrayfun(@(x) reshape([mask{x}([mask{x}.Area] == area(x)).Centroid], 2, []), 1:length(mask), 'uniformoutput',false);
obj.topo.centroid = cell2mat(centroids);

obj.topo.distances = distance(obj.topo.centroid);

if ~isempty(obj.analysis)
    [~, stack] = max( obj.analysis.stack );
    obj.topo.loc = arrayfun(@(x) stack(x) .* (obj.topo.maskNeurons == x), 1:size(obj.twop.deconv,2), 'uniformoutput',false);
    obj.topo.loc = sum( cell2mat( permute(obj.topo.loc, [1 4 3 2]) ), 4 );
end

    function d = distance(coor)
        %distances between centroids
        
        coor = coor .* obj.topo.FOV ./ [size(obj.topo.mimg, 1); size(obj.topo.mimg, 2)];
        d = coor - permute(coor, [1 3 2]);
        d = sqrt( sum(d.^2, 1) );
        d = squeeze(d);
        
    end

end