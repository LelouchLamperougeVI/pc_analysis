%% classify them ensembles... again...
% (run this first)

clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'hiepi_struct', 'pc_list', 'file', 'whole_stack');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'clusts', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'hiepi_struct', 'pc_list', 'file', 'whole_stack');

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
    else%if any(l1{ii} > l_thres)
        iscue1(ii) = 2;
    end
end
iscue2 = zeros(length(l2), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l2) %classify cue/traj ensembles
    if isempty(l2{ii}); continue; end
    temp = any(s2{ii} < cue_centres & e2{ii} > cue_centres & l2{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue2(ii) = 1;
    else%if any(l2{ii} > l_thres)
        iscue2(ii) = 2;
    end
end

% clearvars -except rsc ee iscue1 iscue2


%% build master stacks
% (run this as well)
rest = 2;

hiepi_psth = cat(1, rsc.hiepi_psth{rest}, ee.hiepi_psth{rest});
session = 1:length(hiepi_psth);
session = repelem(session, cellfun(@length, hiepi_psth));
session = session(:);

file = cat(1, rsc.file, ee.file);
ts = cell(length(file), 1);
for ii = 1:length(file)
    temp = load(fullfile(file{ii}, 'analysis.mat'));
    ts{ii} = temp.analysis.behavior.frame_ts;
end
file = repelem(file, cellfun(@length, hiepi_psth));
ts = repelem(ts, cellfun(@length, hiepi_psth));

hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = cellfun(@(x) mean(x, 2), hiepi_psth, 'UniformOutput', false);
hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = fast_smooth(hiepi_psth, 5 / 150 * 50);
raw_hiepi_psth = hiepi_psth;
hiepi_psth = (hiepi_psth - min(hiepi_psth, [], 1)) ./ range(hiepi_psth, 1);

hiepi_z = cat(1, rsc.hiepi_z{rest}, ee.hiepi_z{rest});
hiepi_z = cellfun(@(x) mat2cell(x, size(x, 1), ones(1, size(x, 2))), hiepi_z, 'UniformOutput', false);
hiepi_z = cat(2, hiepi_z{:});
hiepi_z = hiepi_z(:);
% hiepi_z = cellfun(@(x) fast_smooth(x, 4), hiepi_z, 'UniformOutput', false);
hiepi_z = cellfun(@(x) zscore(x(~isnan(x))), hiepi_z, 'UniformOutput', false);

hiepi_struct = cat(1, rsc.hiepi_struct{rest}, ee.hiepi_struct{rest});
hiepi_struct = cat(1, hiepi_struct{:});

whole_stack = cat(2, rsc.whole_stack, ee.whole_stack);
clusts = cat(2, rsc.clusts{rest}, ee.clusts{rest});
clust_stacks_whole = [];
for ii = 1:length(clusts)
    for jj = 1:length(clusts{ii})
        clust_stacks_whole = cat(1, clust_stacks_whole, {whole_stack{ii}(:, clusts{ii}{jj})});
    end
end


%% duration of reactivation events
dur = cellfun(@mean, {hiepi_struct.dur});
figure
boxplot(dur, iscue2)
ylim([0, .6])
kruskalwallis(dur, iscue2)


%%
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'hiepi_struct', 'pc_list', 'file', 'swr_struct');

% classify ensembles
clust_stacks = ee.clust_stacks;

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
    else%if any(l1{ii} > l_thres)
        iscue1(ii) = 2;
    end
end
iscue2 = zeros(length(l2), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l2) %classify cue/traj ensembles
    if isempty(l2{ii}); continue; end
    temp = any(s2{ii} < cue_centres & e2{ii} > cue_centres & l2{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue2(ii) = 1;
    else%if any(l2{ii} > l_thres)
        iscue2(ii) = 2;
    end
end

clearvars -except ee iscue1 iscue2

% params
wdw = 1;
bandwidth = .1;

rest = 2;
swr_on = {ee.swr_struct{rest}(:).swr_on};
swr_on = swr_on(:);

idx = cellfun(@isempty, ee.clusts{rest});
hiepi_on = cell(length(idx), 1);
hiepi_on(~idx) = ee.hiepi_struct{rest}(:);

xcf = [];
for ii = 1:length(swr_on)
    for jj = 1:length(hiepi_on{ii})
        temp = hiepi_on{ii}(jj).on(:)' - swr_on{ii}(:);
        xcf = cat(1, xcf, {temp(:)});
    end
end

xi = linspace(-wdw, wdw, 1e3);
f = cellfun(@(x) ksdensity(x, xi, 'bandwidth', bandwidth), xcf, 'UniformOutput', false);
f = cell2mat(f);
[~, delay] = max(f, [], 2);
delay = xi(delay);

g = repelem(1:2, [length(xi), length(xi)]);
x = repmat(xi(:)', [1, 2]);
y = cat(2, mean(f((iscue2 == 1) & (cellfun(@length, xcf) >= 200) , :)), mean(f((iscue2 == 2) & (cellfun(@length, xcf) >= 200) , :)));
ci = cat(2, sem(f((iscue2 == 1) & (cellfun(@length, xcf) >= 200) , :)), sem(f((iscue2 == 2) & (cellfun(@length, xcf) >= 200) , :)));
g = gramm('x', x, 'y', y, 'ymin', y - ci, 'ymax', y + ci, 'color', g);
g.geom_line;
g.geom_interval;
g.axe_property('ylim', [1e-3, 1.7e-3]);
figure
g.draw;

figure
boxplot(delay(~~iscue2 & (cellfun(@length, xcf) >= 100)), iscue2(~~iscue2 & (cellfun(@length, xcf) >= 100)), 'Orientation', 'horizontal')
xlim([-1, 1]);
p = ranksum(delay((iscue2 == 1) & (cellfun(@length, xcf) >= 100)), delay((iscue2 == 2) & (cellfun(@length, xcf) >= 100)), 'tail', 'left');
title(['ranksum p = ' num2str(p)])

figure
cdfplot(delay(iscue2 == 1))
hold on
cdfplot(delay(iscue2 == 2))
[~, p] = kstest2(delay((iscue2 == 1) & (cellfun(@length, xcf) >= 100)), delay((iscue2 == 2) & (cellfun(@length, xcf) >= 100)), 'Tail', 'larger');
title(['kstest2 p = ' num2str(p)])

groups = cellfun(@length, xcf) >= 100;
subratio = delay(groups);
groups = iscue2(groups);
subratio = subratio(~~groups);
groups = groups(~~groups);
mu = diff(accumarray(groups, subratio, [2, 1], @mean));

p = zeros(length(groups), 1);
for ii = 1:1e5
    idx = randperm(length(groups));
    p(ii) = diff(accumarray(groups(idx), subratio, [2, 1], @mean));
end

figure
histogram(p)
xline(mu)
p = 1 - sum(mu>p)/length(p);
title(['permutation p = ' num2str(p)])