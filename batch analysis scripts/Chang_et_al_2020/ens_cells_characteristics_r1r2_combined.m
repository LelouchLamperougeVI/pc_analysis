clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat');

% classify ensembles
clust_stacks{1} = cat(1, rsc.clust_stacks{1}, ee.clust_stacks{1});
clust_stacks{2} = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});

fr_thres = .5;
traj_thres = 3; %min number of pc per ensemble
l_thres = 30; % length to be considered cue ensemble

l1 = cell(length(clust_stacks{1}), 1); s1=l1; e1=l1;
for c = 1:length(clust_stacks{1})
    stack = clust_stacks{1}{c};
    traj = any(stack > fr_thres, 2);
    [~, starts, ends] = traj_length(traj, 1);

    stack = repmat(stack, 2, 1);
    idx = false(length(starts{1}), 1);
    for t = 1:length(starts{1})
        temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
        temp = any(temp > fr_thres, 1);
        idx(t) = sum(temp) < traj_thres;
    end

    [l1{c}, s1{c}, e1{c}] = traj_length(traj);
    l1{c} = l1{c}{1}(~idx);
    s1{c} = s1{c}{1}(~idx);
    e1{c} = e1{c}{1}(~idx);
end

l2 = cell(length(clust_stacks{2}), 1); s2=l2; e2=l2;
for c = 1:length(clust_stacks{2})
    stack = clust_stacks{2}{c};
    traj = any(stack > fr_thres, 2);
    [~, starts, ends] = traj_length(traj, 1);

    stack = repmat(stack, 2, 1);
    idx = false(length(starts{1}), 1);
    for t = 1:length(starts{1})
        temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
        temp = any(temp > fr_thres, 1);
        idx(t) = sum(temp) < traj_thres;
    end

    [l2{c}, s2{c}, e2{c}] = traj_length(traj);
    l2{c} = l2{c}{1}(~idx);
    s2{c} = s2{c}{1}(~idx);
    e2{c} = e2{c}{1}(~idx);
end

belt;
iscue1 = zeros(length(l1), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l1) %classify cue/traj ensembles
    if isempty(l1{ii}); continue; end
    temp = any(s1{ii} < cue_centres & e1{ii} > cue_centres & l1{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue1(ii) = 1;
    else
        iscue1(ii) = 2;
    end
end
iscue2 = zeros(length(l2), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l2) %classify cue/traj ensembles
    if isempty(l2{ii}); continue; end
    temp = any(s2{ii} < cue_centres & e2{ii} > cue_centres & l2{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue2(ii) = 1;
    else
        iscue2(ii) = 2;
    end
end

clearvars -except rsc ee iscue1 iscue2

width = cat(2, rsc.pc_width, ee.pc_width);
pc_list = cat(2, rsc.pc_list, ee.pc_list);
clusts1 = cat(2, rsc.clusts{1}, ee.clusts{1});
clusts2 = cat(2, rsc.clusts{2}, ee.clusts{2});

num_pf = [];
loc = {};
wid = [];
frac = [];
num_pc_ens = [];
for ii = 1:length(clusts1)
    for jj = 1:length(clusts1{ii})
        frac = cat(1, frac, length(intersect(clusts1{ii}{jj}, pc_list{ii})) ./ length(clusts1{ii}{jj}));
        num_pc_ens = cat(1, num_pc_ens, length(intersect(clusts1{ii}{jj}, pc_list{ii})));
        temp = width{ii}(intersect(clusts1{ii}{jj}, pc_list{ii}));
        num_pf = cat(1, num_pf, mean(cellfun(@(x) size(x, 1), temp)));
        wid = cat(1, wid, mean(cellfun(@(x) mean(x(:, 1)), temp)));
        temp = cat(1, temp{:});
        if isempty(temp)
            loc = cat(1, loc, {[]});
        else
            loc = cat(1, loc, {temp(:, 2)});
        end
    end
end

for ii = 1:length(clusts2)
    for jj = 1:length(clusts2{ii})
        frac = cat(1, frac, length(intersect(clusts2{ii}{jj}, pc_list{ii})) ./ length(clusts2{ii}{jj}));
        num_pc_ens = cat(1, num_pc_ens, length(intersect(clusts2{ii}{jj}, pc_list{ii})));
        temp = width{ii}(intersect(clusts2{ii}{jj}, pc_list{ii}));
        num_pf = cat(1, num_pf, mean(cellfun(@(x) size(x, 1), temp)));
        wid = cat(1, wid, mean(cellfun(@(x) mean(x(:, 1)), temp)));
        temp = cat(1, temp{:});
        if isempty(temp)
            loc = cat(1, loc, {[]});
        else
            loc = cat(1, loc, {temp(:, 2)});
        end
    end
end

clust_stacks = cat(1, rsc.clust_stacks{1}, ee.clust_stacks{1});
r1 = nan(length(clust_stacks), 1);
for ii = 1:length(clust_stacks)
    if isempty(clust_stacks{ii})
        continue
    end
    temp = corr(clust_stacks{ii});
    idx = triu(ones(length(temp)), 1);
    r1(ii) = mean(temp(~~idx));
end
clust_stacks = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});
r2 = nan(length(clust_stacks), 1);
for ii = 1:length(clust_stacks)
    if isempty(clust_stacks{ii})
        continue
    end
    temp = corr(clust_stacks{ii});
    idx = triu(ones(length(temp)), 1);
    r2(ii) = mean(temp(~~idx));
end
r = cat(1, r1, r2);

iscue = cat(1, iscue1, iscue2);

% panel a - mean num place fields
g = iscue(iscue ~= 0);
figure
subplot(1, 2, 1)
cdfplot(num_pf(iscue == 1));
hold on
cdfplot(num_pf(iscue == 2));
[~, p] = kstest2(num_pf(iscue == 1), num_pf(iscue == 2));
title(['kstest p = ' num2str(p)]);

subplot(1, 2, 2)
boxplot(num_pf(iscue ~= 0), g);
ylim([1 2.5])
p = ranksum(num_pf(iscue == 1), num_pf(iscue == 2));
title(['ranksum p = ' num2str(p)]);

% panel b - mean place field width
g = iscue(iscue ~= 0);
figure
subplot(1, 2, 1)
cdfplot(wid(iscue == 1));
hold on
cdfplot(wid(iscue == 2));
[~, p] = kstest2(wid(iscue == 1), wid(iscue == 2));
title(['kstest p = ' num2str(p)]);

subplot(1, 2, 2)
boxplot(wid(iscue ~= 0), g);
ylim([5 35])
p = ranksum(wid(iscue == 1), wid(iscue == 2));
title(['ranksum p = ' num2str(p)]);

% panel c - fraction place cells
figure
boxplot(frac(num_pc_ens > 2), iscue(num_pc_ens > 2))
ylim([0 1])
[p, ~, stats] = kruskalwallis(frac(num_pc_ens > 2), iscue(num_pc_ens > 2))
multcompare(stats, 'ctype', 'bonferroni')

% panel d - place field locations
loc_cue = cat(1, loc{iscue == 1});
loc_traj = cat(1, loc{iscue == 2});

ii = 1e4;
edges = 1:50;
boots = zeros((length(edges) - 1) * 2, ii);
for ii = 1:ii
    boots(:, ii) = cat(1, histcounts(randsample(loc_cue, length(loc_cue), true), edges, 'normalization', 'probability')', histcounts(randsample(loc_traj, length(loc_traj), true), edges, 'normalization', 'probability')');
end
y = cat(1, histcounts(loc_cue, edges, 'normalization', 'probability')', histcounts(loc_traj, edges, 'normalization', 'probability')');
ci = prctile(boots, [2.5, 97.5], 2);
x = repmat((1:(length(edges) - 1))', [2, 1]);
g = repelem((1:2)', [length(edges) - 1; length(edges) - 1]);

g = gramm('x', x, 'y', y, 'ymin', ci(:, 1), 'ymax', ci(:, 2), 'color', g);
g.geom_line;
g.geom_interval;
g.axe_property('ylim', [0 .15]);
figure
g.draw;

% panel e - temporal correlation
g = iscue(iscue ~= 0);
figure
subplot(1, 2, 1)
cdfplot(r(iscue == 1));
hold on
cdfplot(r(iscue == 2));
xlim([-.2 1])
[~, p] = kstest2(r(iscue == 1), r(iscue == 2));
title(['kstest p = ' num2str(p)]);

subplot(1, 2, 2)
boxplot(r(iscue ~= 0), g);
ylim([-.2 1])
p = ranksum(r(iscue == 1), r(iscue == 2));
title(['ranksum p = ' num2str(p)]);


figure
subplot(1, 2, 1)
pie(histcounts(iscue1), [0 1 1])
subplot(1, 2, 2)
pie(histcounts(iscue2), [0 1 1])