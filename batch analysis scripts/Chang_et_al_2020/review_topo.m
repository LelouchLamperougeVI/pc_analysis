clear all

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
% animals = dir(fullfile(root, 'RSC*'));

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

count = 1;

sil = [];
p = [];
d = [];
nni = [];
p_nni = [];
clust = [];

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    disp(['running ' num2str(a) '/' num2str(length(animals))]);
    
    for s = 1:length(sessions)
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.set_ops('e_size',5);
        rest1.set_ops('clust_method','thres');
        rest1.set_ops('sig', .2);
        rest1.remove_mvt;
        rest1.cluster;
        rest1.topography;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
        rest2.topography;
        
        sil = cat(1, sil, {rest1.topo.clust.silhouette, rest2.topo.clust.silhouette});
        p = cat(1, p, {rest1.topo.clust.p, rest2.topo.clust.p});
        d = cat(1, d, {rest1.topo.clust.d, rest2.topo.clust.d});
        nni = cat(1, nni, {rest1.topo.clust.nni, rest2.topo.clust.nni});
        p_nni = cat(1, p_nni, {rest1.topo.clust.p_nni, rest2.topo.clust.p_nni});
        
        topo(count, 1) = rest1.topo;
        topo(count, 2) = rest2.topo;
        
        clust = cat(1, clust, {rest1.ensembles.clust, rest2.ensembles.clust});
        
        count = count + 1;
    end
end

%%
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');

% classify ensembles
clust_stacks{1} = cat(1, rsc.clust_stacks{1}, ee.clust_stacks{1});
clust_stacks{2} = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});
% clust_stacks{1} = rsc.clust_stacks{1};
% clust_stacks{2} = rsc.clust_stacks{2};
clust_stacks{1} = ee.clust_stacks{1};
clust_stacks{2} = ee.clust_stacks{2};

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

clearvars -except iscue1 iscue2


%%
% results1 = load('/mnt/storage/HaoRan/RRR_motor/M2/review_topo.mat');
% results2 = load('/mnt/storage/rrr_magnum/M2/review_topo.mat');
load('/mnt/storage/rrr_magnum/M2/review_topo.mat');

% sil = cat(1, results1.sil, results2.sil);
sil1 = cat(1, sil{:, 1});
sil2 = cat(1, sil{:, 2});
% sil1 = cellfun(@(x) x(:), sil(:, 1), 'UniformOutput', false);
% sil1 = cat(1, sil1{:});
% sil2 = cellfun(@(x) x(:), sil(:, 2), 'UniformOutput', false);
% sil2 = cat(1, sil2{:});
sil = cat(1, sil1, sil2);
sil = cellfun(@mean, sil);

% p_nni = cat(1, results1.p_nni, results2.p_nni);
p_nni = cat(1, p_nni{:, 1}, p_nni{:, 2});
% nni = cat(1, results1.nni, results2.nni);
nni = cat(1, nni{:, 1}, nni{:, 2});

% p = cat(1, results1.p, results2.p);
p = cat(1, p{:, 1}, p{:, 2});

iscue = cat(1, iscue1, iscue2);
% iscue((length(iscue2) + 1):end) = iscue((length(iscue2) + 1):end) + 3;
% [~, p] = kstest2(p(iscue == 1), p(iscue == 2))

% sil = cellfun(@mean, sil2);
% [~, p] = kstest2(sil(iscue2 == 1), sil(iscue2 == 2))
% [~, p] = kstest2(p(iscue2 == 1), p(iscue2 == 2))

figure
subplot(2, 2, 1)
cdfplot(sil(iscue == 1))
hold on
cdfplot(sil(iscue == 2))
xlim([-.5 1])
subplot(2, 2, 2)
cdfplot(nni(iscue == 1))
hold on
cdfplot(nni(iscue == 2))
xlim([0 2])
subplot(2, 2, 3)
boxplot(sil, iscue)
ylim([-.5 1])
subplot(2, 2, 4)
boxplot(nni, iscue)
ylim([0 2])

[~, p] = kstest2(sil(iscue == 1), sil(iscue == 2))
p = ranksum(sil(iscue == 1), sil(iscue == 2))
p = signrank(sil(~~iscue), 0)

[~, p] = kstest2(nni(iscue == 1), nni(iscue == 2))
p = ranksum(nni(iscue == 1), nni(iscue == 2))
p = signrank(nni(~~iscue), 1)

%%
clear all

roots = {'/mnt/storage/HaoRan/RRR_motor/M2/RSC*', '/mnt/storage/rrr_magnum/M2/E*'};
animals = cellfun(@dir, roots, 'UniformOutput', false);
animals = cat(1, animals{:});
animals(~cat(1, animals(:).isdir)) = [];
animals = arrayfun(@(x) fullfile(x.folder, x.name), animals, 'UniformOutput', false);

count = 1;

for a = 1:length(animals)
    sessions = dir(animals{a});
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );

    avg_stack = cell(length(sessions), 1);
    ss = cell(length(sessions), 2);
    for s = 1:length(sessions)
        try
            temp = load(fullfile(animals{a}, sessions{s}, '2', 'Plane1', 'mean_img.mat'));
        catch
            temp = load(fullfile(animals{a}, sessions{s}, '2', 'plane1', 'mean_img.mat'));
        end
        avg_stack{s} = temp.mimg;
        ss{s, 1} = fullfile(animals{a}, sessions{s}, [sessions{s} '_1.abf']);
        ss{s, 2} = fullfile(animals{a}, sessions{s}, [sessions{s} '_3.abf']);
    end
    
    tbl(a).folder = animals{a};
    tbl(a).avg_stack = avg_stack;
    tbl(a).loader = ss;
end


%% window assignment
for ii = 1:length(tbl)
    stack = tbl(ii).avg_stack;
    stack = cat(3, stack{:});
    
    k = 4;
    figure
    for jj = 1:size(stack, 3)
        subplot(k, k, jj);
        imagesc(stack(:, :, jj));
        axis square
        colormap gray
    end
end



%%
load('/mnt/storage/rrr_magnum/M2/review_topo.mat');

centroids = arrayfun(@(x) x.centroid', topo, 'UniformOutput', false);

coor = cell(size(clust));
for ii = 1:size(clust, 1)
    for jj = 1:size(clust, 2)
        coor{ii, jj} = cell(length(clust{ii, jj}), 1);
        for kk = 1:length(clust{ii, jj})
            coor{ii, jj}{kk} = centroids{ii, jj}(clust{ii, jj}{kk}, :);
        end
    end
end

coor1 = cat(1, coor{:, 1});
coor2 = cat(1, coor{:, 2});

nbins = 400;
fov = topo(1, 1).FOV(1);
xy = linspace(0, fov, nbins);
[x, y] = meshgrid(xy, xy);
densities = zeros(nbins, nbins, range(iscue2) + 2);

figure
n = [];
for ii = unique(iscue2)'
    coor = cat(1, coor1{iscue1 == ii}, coor2{iscue2 == ii});
    
    f = kde(coor, [x(:), y(:)], 50, 1, [0, 0; fov, fov]);
    f = reshape(f, nbins, nbins);
    densities(:, :, ii+1) = f;
    n = cat(1, n, size(coor, 1));
end

coor = cat(1, centroids{:, 2});

f = kde(coor, [x(:), y(:)], 50, 1, [0, 0; fov, fov]);
f = reshape(f, nbins, nbins);
densities(:, :, ii+2) = f;
n = cat(1, n, size(coor, 1));

densities = densities ./ sum(densities, [1, 2]);
for ii = 1:size(densities, 3)
    subplot(1, size(densities, 3), ii);
    imagesc('XData', xy, 'YData', xy, 'CData', densities(:, :, ii));
    axis square
    colormap viridis
    title(['n = ' num2str(n(ii))]);
end

temp = log(densities(:, :, 1:end-1) ./ densities(:, :, end));
kl = sum( densities(:, :, 1:end-1) .* temp, [1, 2], 'omitnan');
figure
bar(squeeze(kl));
ylim([0 .1])

%% panels A and B
clear all
load('/mnt/storage/rrr_magnum/M2/review_topo.mat', 'topo', 'clust');
load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'pc_list', 'whole_stack');
load('/mnt/storage/rrr_magnum/M2/ccpc_ll.mat', 'md2');

% reject_idx = [true false true]; % rejected sessions in RSC animals
% reject_idx = repelem(reject_idx, [6 5 1]);
% reject_idx = find(reject_idx);

centroids = arrayfun(@(x) x.centroid', topo, 'UniformOutput', false);
centroids = centroids(:, 2);

stack = [];
coor = [];
ll = [];
for ii = 1:length(pc_list)
    stack = cat(2, stack, whole_stack{ii}(:, pc_list{ii}));
    coor = cat(1, coor, centroids{ii}(pc_list{ii}, :));
    ll = cat(1, ll, md2{ii}.ll(pc_list{ii}, :));
end

fov = topo(1, 1).FOV(1);

nbins = 400;
xy = linspace(0, fov, nbins);
[x, y] = meshgrid(xy, xy);
ratio = -2 .* (ll(:, 2) - ll(:, 1));
f = kde(coor, [x(:), y(:)], 50, ratio, [0, 0; fov, fov]);
f = reshape(f, nbins, nbins);

cm = brewermap([], 'PRGn');

targets = [500, 220
           100, 400
           525, 650
           320, 145];

figure
imagesc('xdata', xy, 'ydata', xy, 'cdata', f);
caxis(max(abs(f(:))) .* [-1 1]);
axis square
colormap(cm);
colorbar
for ii = 1:length(targets)
    rectangle('Position', cat(2, targets(ii, end:-1:1), [100 100]))
end

figure
for ii = 1:length(targets)
    subplot(2, length(targets), ii);
    temp = stack(:, all(coor >= targets(ii, :) & coor <= targets(ii, :) + 100, 2));
    [~, order] = max(temp);
    [~, order] = sort(order);
    imagesc(temp(:, order)');
    colormap jet
    subplot(2, length(targets), ii + length(targets));
    imagesc(corr(temp'));
    axis square
end



%%
clear all
load('/mnt/storage/rrr_magnum/M2/review_topo.mat');
load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'whole_stack');

n = 40;
k = 5;
figure
for ii = 1:length(clust{n, 2})
    temp = whole_stack{n}(:, clust{n, 2}{ii});
    [~, order] = max(temp);
    [~, order] = sort(order);
    
    subplot(k, k, ii)
    imagesc(-temp(:, order)')
    colormap gray
end

figure
for ii = 1:length(clust{n, 2})
    subplot(k, k, ii)
    imagesc(~ismember(topo(n, 2).maskNeurons, clust{n, 2}{ii}));
    colormap gray
    axis square
end

% targets = [4 6 9 10 12];
targets = [4 6];

figure
for ii = 1:length(targets)
    subplot(1, 2, ii)
    imagesc(~ismember(topo(n, 2).maskNeurons, clust{n, 2}{targets(ii)}));
    colormap gray
    axis square
end

figure
subplot(1, 2, 1)
imagesc(topo(n, 2).mimg);
axis square
colormap gray
subplot(1, 2, 2)
imagesc(~topo(n, 2).maskNeurons);
axis square
colormap gray







