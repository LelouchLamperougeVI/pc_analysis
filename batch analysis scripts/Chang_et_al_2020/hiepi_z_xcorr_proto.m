clear all
% load('/mnt/storage/rrr_magnum/M2/swr_stack_2s.mat');
% load('/mnt/storage/rrr_magnum/M2/hiepi.mat');
load('/mnt/storage/rrr_magnum/M2/hiepi4.mat');

fr_thres = .5;
traj_thres = 3; %min number of pc per ensemble
l_thres = 30; % length to be considered cue ensemble

l1 = cell(length(swr_clust_stack{1}), 1); s1=l1; e1=l1;
for s = 1:length(swr_clust_stack{1})
    l1{s} = cell(length(swr_clust_stack{1}{s}), 1);
    for c = 1:length(swr_clust_stack{1}{s})
        stack = swr_clust_stack{1}{s}{c};
        traj = any(stack > fr_thres, 2);
        [~, starts, ends] = traj_length(traj, 1);

        stack = repmat(stack, 2, 1);
        idx = false(length(starts{1}), 1);
        for t = 1:length(starts{1})
            temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
            temp = any(temp > fr_thres, 1);
            idx(t) = sum(temp) < traj_thres;
        end

        [l1{s}{c}, s1{s}{c}, e1{s}{c}] = traj_length(traj);
        l1{s}{c} = l1{s}{c}{1}(~idx);
        s1{s}{c} = s1{s}{c}{1}(~idx);
        e1{s}{c} = e1{s}{c}{1}(~idx);
    end
end

l2 = cell(length(swr_clust_stack{2}), 1); s2=l2; e2=l2;
for s = 1:length(swr_clust_stack{2})
    l2{s} = cell(length(swr_clust_stack{2}{s}), 1);
    for c = 1:length(swr_clust_stack{2}{s})
        stack = swr_clust_stack{2}{s}{c};
        traj = any(stack > fr_thres, 2);
        [~, starts, ends] = traj_length(traj, 1);

        stack = repmat(stack, 2, 1);
        idx = false(length(starts{1}), 1);
        for t = 1:length(starts{1})
            temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
            temp = any(temp > fr_thres, 1);
            idx(t) = sum(temp) < traj_thres;
        end

        [l2{s}{c}, s2{s}{c}, e2{s}{c}] = traj_length(traj);
        l2{s}{c} = l2{s}{c}{1}(~idx);
        s2{s}{c} = s2{s}{c}{1}(~idx);
        e2{s}{c} = e2{s}{c}{1}(~idx);
    end
end

belt;
iscue1 = cell(length(l1), 1);
for ii = 1:length(l1) %classify cue/traj ensembles
    iscue1{ii} = zeros(length(l1{ii}), 1);
    for jj = 1:length(l1{ii})
        if isempty(l1{ii}{jj}); continue; end
        temp = any(s1{ii}{jj} < cue_centres & e1{ii}{jj} > cue_centres & l1{ii}{jj} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
        if temp
            iscue1{ii}(jj) = 1;
        else
            iscue1{ii}(jj) = 2;
        end
    end
end
iscue2 = cell(length(l2), 1);
for ii = 1:length(l2) %classify cue/traj ensembles
    iscue2{ii} = zeros(length(l2{ii}), 1);
    for jj = 1:length(l2{ii})
        if isempty(l2{ii}{jj}); continue; end
        temp = any(s2{ii}{jj} < cue_centres & e2{ii}{jj} > cue_centres & l2{ii}{jj} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
        if temp
            iscue2{ii}(jj) = 1;
        else
            iscue2{ii}(jj) = 2;
        end
    end
end