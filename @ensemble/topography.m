function topography(obj)
% Tales from Topographic Oceans :)
% Obscure reference from the 70's...

disp('Running superclass topography method.');
topography@LFP(obj);

if isempty(obj.ensembles.clust)
    warning('Data either not clustered, or no ensemble detected. Skipping subclass subroutine.');
    obj.topo.clust.silhouette = {};
    obj.topo.clust.p = [];
    obj.topo.clust.d = {};
    obj.topo.clust.nni = [];
    obj.topo.clust.p_nni = [];
    return
end

% obj.topo.clust.vertices = cellfun(@vertices, obj.topo.clust.centroid, 'uniformoutput',false);

% % get silhouette score between clusters
% X = obj.topo.centroid(:, cell2mat(obj.ensembles.clust))';
% k = arrayfun(@(x) x .* ones(1, length(obj.ensembles.clust{x})), 1:length(obj.ensembles.clust), 'uniformoutput', false);
% k = cell2mat(k)';
% s = silhouette(X, k, 'sqEuclidean');
% obj.topo.clust.silhouette = arrayfun(@(x) s(k==x), unique(k), 'uniformoutput', false );
% % silhouette over z-aspect as well
% X = [obj.topo.centroid(:, cell2mat(obj.ensembles.clust))', obj.twop.planes.depth(cell2mat(obj.ensembles.clust))'];
% s = silhouette(X, k, 'sqEuclidean');
% obj.topo.clust.zsilhouette = arrayfun(@(x) s(k==x), unique(k), 'uniformoutput', false );
%
% obj.topo.clust.k = k;

d = obj.topo.centroid';
k = zeros(size(d, 1), 1);
for ii = 1:length(obj.ensembles.clust)
    k(obj.ensembles.clust{ii}) = ii;
end
md = silhouette(d, k);
obj.topo.clust = md;

function v = vertices(cent)
% C x P x V matrix
% C: xy coor; P: coor pair; V: vertices
idx = nchoosek(1:size(cent,2), 2);
v = zeros(2, 2, size(idx,1));
for ii = 1:size(idx,1)
    v(:,:,ii) = cent(:,idx(ii,:));
end

function md = silhouette(d, k)
% Silhouette scores for topography analysis
% Unlike the actual silhouette score, where the minimum mean distance is
% taken between the target cluster and all other clusters, the other
% clusters in question is considered to be all neurons not part of the
% current ensemble.
%
% Inputs
%   d       m-by-n matrix of coordinates with rows as observations and
%           columns as dimensions
%   k       m-by-1 vector of cluster labels. 0 are non-ensemble neurons
% Outputs
%   s       silhouette scores for each neuron in each ensemble

it = 1e3;

d = permute(d, [1, 3, 2]) - permute(d, [3, 1, 2]);
d = sqrt(sum(d .^ 2, 3)); % Euclidean
% d = sum(d .^ 2, 3); % squared Euclidean

s = cell(range(k), 1);
p = zeros(range(k), 1);
dist = cell(range(k), 1);
nni = zeros(range(k), 1); % nearest neighbour index
p_nni = zeros(range(k), 1);
for ii = 1:range(k)
    idx = k == ii;
    a = sum(d(idx, idx, :)) ./ (sum(idx) - 1);
    b = sum(d(~idx, idx, :)) ./ sum(~idx);
    temp = (b - a) ./ max(cat(1, a, b));
    s{ii} = temp(:);
    
    mu = zeros(it, 1);
    mu_d = zeros(it, 1);
    for jj = 1:it
        sample = randsample(length(d), sum(idx));
        d_sample = d(sample, sample);
        d_sample(~~eye(length(d_sample))) = inf;
        mu_d_sample = min(d_sample);
        mu_d(jj) = mean(mu_d_sample);
        temp = d_sample(triu(true(length(d_sample)), 1));
        mu(jj) = median(temp);
    end
    d_sample = d(idx, idx);
    d_sample(~~eye(length(d_sample))) = inf;
    mu_d_sample = min(d_sample);
    nni(ii) = mean(mu_d_sample) / mean(mu_d);
    p_nni(ii) = sum( mean(mu_d_sample) >= mu_d ) ./ it;
    temp = d_sample(triu(true(length(d_sample)), 1));
    dist{ii} = temp(:);
    temp = median(temp);
    p(ii) = sum(temp >= mu) ./ it;
end

md.silhouette = s;
md.p = p;
md.d = dist;
md.nni = nni;
md.p_nni = p_nni;