function topography(obj)
% Tales from Topographic Oceans :)
% Obscure reference from the 70's; don't worry if you don't get it

if isempty(obj.clust)
    error('data was not clustered');
end

topography@LFP(obj, obj.ops.FOV);

obj.topo.clust.masks = cellfun(@(x) ismember(obj.topo.maskNeurons, x), obj.clust, 'uniformoutput',false);

obj.topo.clust.centroid = cellfun(@regionprops, obj.topo.clust.masks, 'uniformoutput',false);
obj.topo.clust.centroid = arrayfun(@(x) reshape([obj.topo.clust.centroid{x}.Centroid], 2, []), 1:length(obj.topo.clust.centroid), 'uniformoutput',false);

obj.topo.clust.vertices = cellfun(@vertices, obj.topo.clust.centroid, 'uniformoutput',false);

obj.topo.clust.distances = cellfun(@distance, obj.topo.clust.centroid, 'uniformoutput',false);
obj.topo.clust.mu_d = cellfun(@(x) mean(distance(x)), obj.topo.clust.centroid, 'uniformoutput',false);

obj.topo.clust.masks = permute(1:length(obj.clust), [3 1 2]) .* cell2mat( permute(obj.topo.clust.masks, [3 1 2]) );
obj.topo.clust.masks = sum( obj.topo.clust.masks, 3);


    function d = distance(coor)
        %distances between centroids
        coor = coor .* obj.ops.FOV ./ size(obj.topo.mimg)';
        d = coor - permute(coor, [1 3 2]);
        d = sqrt( sum(d.^2, 1) );
        d = squeeze(d);
        d = tril(d, -1);
        d = d( tril(true(size(d)), -1) );
    end

    function v = vertices(cent)
        % C x P x V matrix
        % C: xy coor; P: coor pair; V: vertices
        idx = nchoosek(1:size(cent,2), 2);
        v = zeros(2, 2, size(idx,1));
        for ii = 1:size(idx,1)
            v(:,:,ii) = cent(:,idx(ii,:));
        end
    end

end