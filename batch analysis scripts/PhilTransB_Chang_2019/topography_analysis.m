%% Analysis of cluster topography for RRR (Phil Trans B)
clear all

path = '/mnt/storage/HaoRan/RRR';
animal = 'RSC036';
date = '2017_09_11';
sess = '3';

%% Single session analysis
% Prototyping with a single session

ass = ensemble(fullfile(path, animal, date, [date '_' sess '.abf']));
ass.remove_mvt;
ass.set_ops('e_size',10);
ass.set_ops('clust_method','thres');
ass.set_ops('order','cluster');
ass.cluster;
ass.topography;

%% Silhouette score
X = cell2mat(ass.topo.clust.centroid)';
k = arrayfun(@(x) x .* ones(1, size(ass.topo.clust.centroid{x}, 2)), 1:length(ass.topo.clust.centroid), 'uniformoutput', false);
k = cell2mat(k)';

% X = ass.topo.centroid';
% k = ones(size(X, 1), 1);
% for ii = 1:length(ass.clust)
%     k(ass.clust{ii}) = ii + 1;
% end

figure
[s, h] = silhouette(X, k, 'sqEuclidean');

%% Distance distributions
d_all = squareform(ass.topo.distances)';
d_clust = cell2mat(ass.topo.clust.distances');
k = arrayfun(@(x) x .* ones(1, length(ass.topo.clust.distances{x})), 1:length(ass.topo.clust.distances), 'uniformoutput', false);
k = cell2mat(k)';
k = [k; zeros(length(d_all), 1)];
d = [d_clust; d_all];

[p, ~, stats] = kruskalwallis(d, k);
multcompare(stats)


%% Batch analyse everything!!!
clear all

path = '/mnt/storage/HaoRan/RRR';
animal = {'RSC036', 'RSC037', 'RSC038'};

sil1 = [];
sil2 = [];
d1 = [];
k1 = [];

for ii = 1:length(animal)
    root = dir(fullfile(path, animal{ii}));
    for jj = 3:length(root)
        try
            ass = ensemble(fullfile(path, animal{ii}, root(jj).name, [root(jj).name '_1.abf']));
            ass.remove_mvt;
            ass.set_ops('e_size',10);
            ass.set_ops('clust_method','thres');
            ass.set_ops('order','cluster');
            ass.cluster;
            ass.topography;

            sil1 = [sil1; cell2mat(ass.topo.clust.silhouette)];
            d_all = squareform(ass.topo.distances)';
            d_clust = cell2mat(ass.topo.clust.distances');
            k = arrayfun(@(x) ones(1, length(ass.topo.clust.distances{x})), 1:length(ass.topo.clust.distances), 'uniformoutput', false);
            k = cell2mat(k)';
            k1 = [k1; k; zeros(length(d_all), 1)];
            d1 = [d1; d_clust; d_all];

            ass = ensemble(fullfile(path, animal{ii}, root(jj).name, [root(jj).name '_3.abf']));
            ass.remove_mvt;
            ass.set_ops('e_size',10);
            ass.set_ops('clust_method','thres');
            ass.set_ops('order','cluster');
            ass.cluster;
            ass.topography;

            sil2 = [sil2; cell2mat(ass.topo.clust.silhouette)];
            d_clust = cell2mat(ass.topo.clust.distances');
            k = arrayfun(@(x) 2 .* ones(1, length(ass.topo.clust.distances{x})), 1:length(ass.topo.clust.distances), 'uniformoutput', false);
            k = cell2mat(k)';
            k1 = [k1; k];
            d1 = [d1; d_clust];
        catch
        end
    end
end

%%
[p, ~, stats] = kruskalwallis(d1, k1);
multcompare(stats)
figure
boxplot(d1, k1, 'plotstyle', 'compact');

figure
cdfplot(sil1)
hold on
cdfplot(sil2)






