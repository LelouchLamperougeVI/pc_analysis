clear all

circ_dist = @(X0, X) min(cat(3, mod(X0 - X, 50), mod(X - X0, 50)), [], 3) .* 150 ./ 50;

thres = 3; %min number of pc per ensemble

loc1 = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'loc', 'loc_clust', 'clusts', 'pc_list');
loc2 = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'loc', 'loc_clust', 'clusts', 'pc_list');
loc = cat(2, loc1.loc, loc2.loc);
loc_clust{1} = cat(1, loc1.loc_clust{1}, loc2.loc_clust{1});
loc_clust{2} = cat(1, loc1.loc_clust{2}, loc2.loc_clust{2});
clusts{1} = cat(2, loc1.clusts{1}, loc2.clusts{1});
clusts{2} = cat(2, loc1.clusts{2}, loc2.clusts{2});
clusts_len{1} = cellfun(@length, clusts{1});
clusts_len{2} = cellfun(@length, clusts{2});
pc_list = cat(2, loc1.pc_list, loc2.pc_list);

d1 = []; d2 = [];
loc_cum1 = []; loc_cum2 = [];

s1 = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = arrayfun(@(x) accumarray(loc_clust{1}{x} + 1, 1, [clusts_len{1}(x) + 1, 1]), 1:length(loc_clust{1}), 'uniformoutput', false);
% g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 2:length(g{ii})
        if g{ii}(jj) < thres
            d1 = cat(1, d1, nan);
            continue
        end
        d = loc{ii}(loc_clust{1}{ii} == (jj - 1));
        loc_cum1 = cat(1, loc_cum1, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d1 = cat(1, d1, mean(d));
    end
end
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
g = arrayfun(@(x) ismember(loc_clust{1}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s1 = arrayfun(@(x) s1{x}(g{x}), 1:length(s1), 'uniformoutput', false);

s2 = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = arrayfun(@(x) accumarray(loc_clust{2}{x} + 1, 1, [clusts_len{2}(x) + 1, 1]), 1:length(loc_clust{2}), 'uniformoutput', false);
% g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 2:length(g{ii})
        if g{ii}(jj) < thres
            d2 = cat(1, d2, nan);
            continue
        end
        d = loc{ii}(loc_clust{2}{ii} == (jj - 1));
        loc_cum2 = cat(1, loc_cum2, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d2 = cat(1, d2, mean(d));
    end
end
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
g = arrayfun(@(x) ismember(loc_clust{2}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s2 = arrayfun(@(x) s2{x}(g{x}), 1:length(s2), 'uniformoutput', false);

figure
cdfplot(cell2mat(s1'));
hold on
cdfplot(cell2mat(s2'));

[~, p] = kstest2(cell2mat(s1'), cell2mat(s2'));
disp(['kstest silhouette: ' num2str(p)])

figure
boxplot([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)])

% p = kruskalwallis([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)]);
p = ranksum(cell2mat(s1'), cell2mat(s2'));
disp(['ranksum silhouette: ' num2str(p)])

figure;
violin({d1, d2}, 'labels', {'rest1', 'rest2'}, 'bandwidth', 2)
% violin({d1, d2}, 'labels', {'rest1', 'rest2'}, 'scatter')

% p = kruskalwallis([d1; d2], [ones(length(d1), 1); 2 .* ones(length(d2), 1)]);
p = ranksum(d1, d2);
disp(['ranksum dist: ' num2str(p)])

figure;
histogram(loc_cum1, 50, 'normalization', 'probability')
hold on
histogram(loc_cum2, 50, 'normalization', 'probability')