%% classify them ensembles... again...

clear all
load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'whole_stack', 'clusts');

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

iscue1 = mat2cell(iscue1, cellfun(@length, clusts{1}));
iscue2 = mat2cell(iscue2, cellfun(@length, clusts{2}));

clearvars -except iscue1 iscue2 whole_stack clusts

%%
root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals([animals.isdir]).name};

anal_f = cell(length(animals), 1); % grosse cochonne, tu pensais a quoi? ;)
daydate = cell(length(animals), 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    anal_f{a} = cellfun(@(x) fullfile(root, animals{a}, x, 'analysis.mat'), sessions, 'uniformoutput', false);
    
    daydate{a} = cellfun(@(x) datetime(x, 'InputFormat', 'yyyy_MM_dd'), sessions);
end

idx = cellfun(@length, daydate) > 2;
for ii = find(idx(:)')
    daydate{ii} = [0 cumsum(split((caldiff(daydate{ii})), {'days'}))];
end
for ii = find(~idx(:)')
    if ~isempty(daydate{ii})
        daydate{ii} = 0;
    else
        daydate{ii} = [];
    end
    anal_f{ii} = '';
end


%% Consecutive days mode (this block is mutually exclusive with next block)
target_days = [0 1];

idx = cell(length(daydate), 1);
for ii = 1:length(daydate)
    idx{ii} = [(1:(length(daydate{ii})-1))', (2:length(daydate{ii}))'];
    anal_f{ii} = anal_f{ii}(idx{ii});
    temp = iscue2(sum(cellfun(@length, daydate(1:(ii-1))))+1 : sum(cellfun(@length, daydate(1:ii)))-1)';
    anal_f{ii}(~cellfun(@any, temp), :) = [];
    idx{ii}(~cellfun(@any, temp), :) = [];
    idx{ii} = idx{ii} + sum(cellfun(@length, daydate(1:(ii-1))));
end

idx = cell2mat(idx);
anal_f = cat(1, anal_f{:});


%% Register them masks and find overlapping ROIs
% Quand la realite frappe a la porte... t'as rien a faire que de te
% preparer a te faire mettre :)

olap_thres = .5;

load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'masks');
masks = masks{2};

jaccard = [];
% this code was stolen from multiplane/register.m
for s = 1:size(idx, 1)
    for r = 1:length(target_days)
        for c = 1:length(target_days)
            maskReg = masks{idx(s, r)};
            maskA = regmasks(maskReg, masks{idx(s, c)});
            
            stack = maskReg == permute(1:max(maskReg(:)), [1 3 2]); % tried to be lazy, ended up overcomplicating things...
            stack = maskA .* stack;
            stack = reshape(stack, [numel(maskA) size(stack,3)]); % le gros machin complique a transformations des matrices irregulieres
            stack = arrayfun(@(x) accumarray(stack(:, x)+1, 1, [max(maskA(:))+1 1]), 1:size(stack,2), 'uniformoutput',false);
            stack = cell2mat(stack);
            stack = stack(2:end,:);
            
            movPix = accumarray(maskReg(:)+1, 1, [max(maskReg(:))+1 1]);
            movPix = movPix(2:end);
            
            fixPix = accumarray(maskA(:)+1, 1, [max(maskA(:))+1 1]); % number of pixels for fixed ROIs
            fixPix = fixPix(2:end);
            
            stack = stack ./ ( fixPix + movPix' - stack );
            jaccard = cat(1, jaccard, stack(stack > 0));
        end
    end
end
sess_idx = idx;


%% modelling jaccard overlap probabilities
% I.e., my bastardized version of Sheintuch 2017

clear all
load('/mnt/storage/rrr_magnum/M2/roi_jaccard.mat');
k = kmeans(jaccard, 2);

if mean(jaccard(k == 1)) > mean(jaccard(k == 2))
    k = ~(k - 1) + 1;
end

pl = betafit(jaccard(k == 1));
ph = betafit(jaccard(k == 2));

x = linspace(0, 1, 1e5);
edges = linspace(0, 1, 1e2);

figure
plot(x, betapdf(x, pl(1), pl(2)))
hold on
plot(x, betapdf(x, ph(1), ph(2)))
histogram(jaccard, edges, 'Normalization', 'pdf');
xlim([0 1])
ylim([0 20])
title(['type I: ' num2str(1 - betacdf(.5, pl(1), pl(2))) ' type II:' num2str(betacdf(.5, ph(1), ph(2)))]);