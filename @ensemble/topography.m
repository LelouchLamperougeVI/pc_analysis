function topography(obj)
% Tales from Topographic Oceans :)
% Obscure reference from the 70's; don't worry if you don't get it

if isempty(obj.clust)
    warning('data either not clustered, or no ensemble detected');
    return
end

topography@LFP(obj, obj.ops.FOV);

obj.topo.clust.masks = cellfun(@(x) ismember(obj.topo.maskNeurons, x), obj.clust, 'uniformoutput',false);

obj.topo.clust.centroid = cellfun(@regionprops, obj.topo.clust.masks, 'uniformoutput',false);
obj.topo.clust.centroid = arrayfun(@(x) reshape([obj.topo.clust.centroid{x}.Centroid], 2, []), 1:length(obj.topo.clust.centroid), 'uniformoutput',false);

obj.topo.clust.vertices = cellfun(@vertices, obj.topo.clust.centroid, 'uniformoutput',false);

obj.topo.clust.distances = cellfun(@distance, obj.topo.clust.centroid, 'uniformoutput',false);
obj.topo.clust.mu_d = cellfun(@(x) mean(distance(x)), obj.topo.clust.centroid);

% the following block was hastily implemented to test the "clusteredness"
% of each ensemble in physical space via a permutation test
shuffles = 1000;
mu_null = zeros(shuffles, length(obj.topo.clust.centroid));
centroid = obj.topo.clust.centroid;
centroids = obj.topo.centroid;
fhandle = @distance;
parfor kk = 1:shuffles
    temp = cellfun(@(x) randperm(size(centroids,2), size(x,2)), centroid, 'uniformoutput',false);
    temp = cellfun(@(x) centroids(:, x), temp, 'uniformoutput',false);
    mu_null(kk, :) = cellfun(@(x) mean( feval( fhandle, x)), temp);
end
obj.topo.clust.mu_null = mu_null;
obj.topo.clust.pval = sum( obj.topo.clust.mu_d > mu_null ) ./ shuffles;
% end block

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