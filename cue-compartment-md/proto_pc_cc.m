clear all

traj_thres = .2;
root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));
% root = '/mnt/storage/rrr_magnum/M2';
% animals = dir(fullfile(root, 'E*'));
animals = {animals.name};

clust_stacks = cell(2,1);
trajectories = cell(2,1);
origin = cell(2,1);
loc_clust = cell(2,1);
clusts = cell(2,1);
pev = cell(2, 1);
file_ids = cell(2, 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.set_ops('e_size',5);
        rest1.set_ops('clust_method','thres');
        rest1.set_ops('sig', .2);
        rest1.remove_mvt;
        rest1.cluster;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
        
        file_ids{1} = cat(1, file_ids{1}, {fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1'])});
        file_ids{2} = cat(1, file_ids{2}, {fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3'])});
        
        % v checkpoint
    
        clusts{1} = cat(2, clusts{1}, {rest1.ensembles.clust});
        clusts{2} = cat(2, clusts{2}, {rest2.ensembles.clust});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        
        analysis = rest1.analysis;
        L = length(analysis.Pi);
        thres=noRun(analysis.behavior.unit_vel);
        thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
        n = analysis.original_deconv;
        dt = .3; % 300 milliseconds
        dt = round(analysis.fs * dt);
        x = analysis.behavior.unit_pos;
            
        for c = 1:length(rest1.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest1.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            clust_stacks{1} = cat(1, clust_stacks{1}, {stack});
            trajectories{1} = cat(2, trajectories{1}, any(stack > traj_thres, 2));
            origin{1} = cat(1, origin{1}, {fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf'])});
            
            try
                md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false);
                pev{1} = cat(1, pev{1}, {[md.pc.pev(:), md.cc.pev(:)]});
            catch
                pev{1} = cat(1, pev{1}, {[]});
            end
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{1} = cat(1, loc_clust{1}, {temp});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        
        analysis = rest2.analysis;
        L = length(analysis.Pi);
        thres=noRun(analysis.behavior.unit_vel);
        thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
        n = analysis.original_deconv;
        dt = .3; % 300 milliseconds
        dt = round(analysis.fs * dt);
        x = analysis.behavior.unit_pos;
        
        for c = 1:length(rest2.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest2.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            clust_stacks{2} = cat(1, clust_stacks{2}, {stack});
            trajectories{2} = cat(2, trajectories{2}, any(stack > traj_thres, 2));
            origin{2} = cat(1, origin{2}, {fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf'])});
            
            try
                md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false);
                pev{2} = cat(1, pev{2}, {[md.pc.pev(:), md.cc.pev(:)]});
            catch
                pev{2} = cat(1, pev{2}, {[]});
            end
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{2} = cat(1, loc_clust{2}, {temp});
    end
end


%%
clear all
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/proto.mat');
ee = load('/mnt/storage/rrr_magnum/M2/proto.mat');
clust_stacks{1} = [rsc.clust_stacks{1}; ee.clust_stacks{1}];
clust_stacks{2} = [rsc.clust_stacks{2}; ee.clust_stacks{2}];
pev{1} = [rsc.pev{1}; ee.pev{1}];
pev{2} = [rsc.pev{2}; ee.pev{2}];

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



%%
% cellfun(@(x) size(x, 2), clust_stacks{2});

n = 0;
idx = cellfun(@(x) size(x, 1) >= n, pev{2});

ev = pev{2}(idx);
ev = cellfun(@sqrt, ev, 'UniformOutput', false);
temp = cellfun(@(x) median(x, 1, 'omitnan'), ev, 'UniformOutput', false);
stack = clust_stacks{2}(idx);
stack(cellfun(@isempty, temp)) = [];
ev(cellfun(@isempty, temp)) = [];
k = iscue2(idx) + 1;
k(cellfun(@isempty, temp)) = [];
temp = cell2mat( temp(~cellfun(@isempty, temp)) );

for ii = 1:length(stack)
    [~, idx] = max(stack{ii}, [], 1);
    [~, idx] = sort(idx);
    stack{ii} = stack{ii}(:, idx);
end

k = kmeans(temp, 3);
% k(:) = 1;
% k(temp(:, 1) > .04) = 2;
% k(temp(:, 1) < .04 & temp(:, 2) > .19) = 3;

figure
hold on
n = 1; scatter(temp(k==n,1), temp(k==n,2))
n = 2; scatter(temp(k==n,1), temp(k==n,2))
n = 3; scatter(temp(k==n,1), temp(k==n,2))

for ii = 1:length(unique(k))
    list = find(k == ii);
    for jj = 1:length(list)
        if ~mod(jj-1, 25)
            figure
        end
        subplot(5, 5, mod(jj-1, 25)+1)
        imagesc(stack{list(jj)}')
        title(num2str(list(jj)))
    end
end


%% Troubleshooting
clear all

rest = ensemble('/mnt/md0/Data/RSC_M2/RSC037/2017_09_13/2017_09_13_3.abf');
% rest = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_09_13/2017_09_13_3.abf');
rest.set_ops('e_size',5);
rest.set_ops('clust_method','thres');
rest.set_ops('sig', .2);
rest.remove_mvt;
rest.cluster;

%%
analysis = rest.analysis;
thres=noRun(analysis.behavior.unit_vel);
thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
n = analysis.original_deconv;
dt = .3; % 300 milliseconds
dt = round(analysis.fs * dt);
x = analysis.behavior.unit_pos;
list = rest.ensembles.clust{13};

% ii = 1;
% parfor ii = 1:100
    md(ii) = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false);
    ii
% end




