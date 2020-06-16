function topography(obj)
% does a bunch of stuff
% the important part is the obj.topo.loc which converts cell number indices
% to positional bin number

distance = @(coor) squeeze( sqrt( sum( (coor - permute(coor, [1 3 2])) .^ 2, 1) ) );

obj.twop.planes.depth = (obj.twop.planes.planes - 1) .* obj.twop.planes.stepsize;
depth = zeros(length(obj.twop.planes.plane_members), 1);
for ii = 1:length(obj.twop.planes.planes)
    depth(obj.twop.planes.plane_members == obj.twop.planes.planes(ii)) = obj.twop.planes.depth(ii);
end
obj.twop.planes.depth = depth';

mask = arrayfun(@(x) regionprops(any(obj.topo.maskNeurons == x, 3)), 1:size(obj.twop.deconv, 2), 'uniformoutput', false);
area = cellfun(@(x) max([x.Area]), mask);
centroids = arrayfun(@(x) reshape([mask{x}([mask{x}.Area] == area(x)).Centroid], 2, []), 1:length(mask), 'uniformoutput',false);
obj.topo.centroid = cell2mat(centroids);
obj.topo.centroid = obj.topo.centroid ./ [size(obj.topo.mimg, 1); size(obj.topo.mimg, 2)] .* obj.topo.FOV;

obj.topo.distances = distance(obj.topo.centroid); % distance over only XY aspect
obj.topo.zdistances = distance(cat(1, obj.topo.centroid, obj.twop.planes.depth));

if ~isempty(obj.analysis)
    [~, stack] = max( obj.analysis.stack );
    obj.topo.loc = arrayfun(@(x) stack(x) .* (obj.topo.maskNeurons == x), 1:size(obj.twop.deconv,2), 'uniformoutput',false);
    obj.topo.loc = sum( cell2mat( permute(obj.topo.loc, [1 4 3 2]) ), 4 );
end