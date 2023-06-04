%% Run main analysis

clear all

root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));

animals = {animals.name};
file = {}; % list of sessions

centroids = [];
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        file = cat(1, file, fullfile(root, animals{a}, sessions{s}));
        
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.topography;
        
        centroids = cat(1, centroids, {rest1.topo.centroid'});
    end
end

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};
file = {}; % list of sessions

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        file = cat(1, file, fullfile(root, animals{a}, sessions{s}));
        
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.topography;
        
        centroids = cat(1, centroids, {rest1.topo.centroid'});
    end
end

%% load up data
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');

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

clearvars -except rsc ee iscue1 iscue2

clusts1 = cat(1, rsc.clusts{1}', ee.clusts{1}');
clusts2 = cat(1, rsc.clusts{2}', ee.clusts{2}');

load('/mnt/storage/rrr_magnum/M2/topography.mat');


%% Calculate variances
v = [];
v_null = [];
for ii = 1:length(clusts1)
    for jj = 1:length(clusts1{ii})
        temp = centroids{ii}(clusts1{ii}{jj}, :);
        v = cat(1, v, var(temp));
        
        temp = centroids{ii};
        v_null = cat(1, v_null, var(temp));
    end
end

for ii = 1:length(clusts2)
    for jj = 1:length(clusts2{ii})
        temp = centroids{ii}(clusts2{ii}{jj}, :);
        v = cat(1, v, var(temp));
        
        temp = centroids{ii};
        v_null = cat(1, v_null, var(temp));
    end
end

fratio = v ./ v_null;
iscue = cat(1, iscue1, iscue2);

figure
g = repmat(1:2, length(iscue), 1);
x = repmat(iscue, 1, 2);
g = gramm('x', x(:), 'y', fratio(:), 'color', g(:));
g.stat_violin('half', true);
g.axe_property('xlim', [-1 3], 'ylim', [0 3.5]);
g.draw;

figure
g = repmat(1:2, length(iscue), 1);
x = repmat(iscue, 1, 2);
g = gramm('x', x(:), 'y', fratio(:), 'color', g(:));
g.stat_boxplot;
g.axe_property('xlim', [-1 3], 'ylim', [0 3.5]);
g.draw;


%% mimg
load('/mnt/storage/rrr_magnum/M2/EE006/2019_06_22/2/plane1/masks_neurons.mat');
a = clusts2{85}{4};
b = clusts2{85}{8};

mimg = zeros(size(maskNeurons, 1), size(maskNeurons, 2), 3);
mimg(:, :, 1) = ~~maskNeurons;
mimg(:, :, 2) = ismember(maskNeurons, a);
mimg(:, :, 3) = ismember(maskNeurons, b);
mimg(:, :, 1) = mimg(:, :, 1) - any(mimg(:, :, 2:3), 3);
mimg = mimg + ~any(mimg, 3);

figure
imshow(mimg)



%% PC topography hypothesis
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'whole_stack', 'pc_list');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'whole_stack', 'pc_list');
pc_list = cat(1, rsc.pc_list', ee.pc_list');
stack = cat(1, rsc.whole_stack', ee.whole_stack');
load('/mnt/storage/rrr_magnum/M2/topography.mat');

pc_d = [];
topo_d = [];
for ii = 1:length(stack)
    temp = stack{ii}(:, pc_list{ii});
    [~, temp] = max(temp, [], 1);
    temp = abs(temp(:) - temp(:)');
    temp = cat(3, temp, 50 - temp);
    temp = min(temp, [], 3);
    idx = triu(true(length(temp)), 1);
    pc_d = cat(1, pc_d, temp(idx));
    
    temp = centroids{ii}(pc_list{ii}, :);
    temp = sqrt( (temp(:, 1) - temp(:, 1)').^2 + (temp(:, 2) - temp(:, 2)').^2 );
    idx = triu(true(length(temp)), 1);
    topo_d = cat(1, topo_d, temp(idx));
end

cm = accumarray([pc_d, discretize(topo_d, linspace(0, 1200, 201))] + 1, 1);
cm = fast_smooth(cm, 5, 2);
imagesc((cm ./ sum(cm, 2))')
colormap jet
















