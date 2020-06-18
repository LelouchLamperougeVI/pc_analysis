function topography(obj)
% Tales from Topographic Oceans :)
% Obscure reference from the 70's...

disp('Running superclass topography method.');
topography@LFP(obj);

if isempty(obj.clust)
    warning('Data either not clustered, or no ensemble detected. Skipping subclass subroutine.');
    return
end

% obj.topo.clust.vertices = cellfun(@vertices, obj.topo.clust.centroid, 'uniformoutput',false);

% get silhouette score between clusters
X = obj.topo.centroid(:, cell2mat(obj.clust))';
k = arrayfun(@(x) x .* ones(1, length(obj.clust{x})), 1:length(obj.clust), 'uniformoutput', false);
k = cell2mat(k)';
s = silhouette(X, k, 'sqEuclidean');
obj.topo.clust.silhouette = arrayfun(@(x) s(k==x), unique(k), 'uniformoutput', false );
% silhouette over z-aspect as well
X = [obj.topo.centroid(:, cell2mat(obj.clust))', obj.twop.planes.depth(cell2mat(obj.clust))'];
s = silhouette(X, k, 'sqEuclidean');
obj.topo.clust.zsilhouette = arrayfun(@(x) s(k==x), unique(k), 'uniformoutput', false );

obj.topo.clust.k = k;

function v = vertices(cent)
% C x P x V matrix
% C: xy coor; P: coor pair; V: vertices
idx = nchoosek(1:size(cent,2), 2);
v = zeros(2, 2, size(idx,1));
for ii = 1:size(idx,1)
    v(:,:,ii) = cent(:,idx(ii,:));
end