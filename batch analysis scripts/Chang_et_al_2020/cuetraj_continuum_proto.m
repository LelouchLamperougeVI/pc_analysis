%% run batch analysis
clear all

traj_thres = .2;

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
% animals = dir(fullfile(root, 'RSC*'));

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

md1 = []; md2 = [];
count = 1;
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
%         rest1.hPICA;
        md1{count} = ccpc_ll(rest1);
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
%         rest2.hPICA;
        md2{count} = ccpc_ll(rest2);
        
        count = count + 1;
    end
end



%% classify ensembles
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

%% quantify likelihood ratios between cue/pc models
ee_md = load('/mnt/storage/rrr_magnum/M2/ccpc_ll.mat'); % this one works
rsc_md = load('/mnt/storage/HaoRan/RRR_motor/M2/ccpc_ll.mat');
% ee_md = load('/mnt/storage/rrr_magnum/M2/ccpc_ll_2.mat'); % with baseline and sd contraint
% rsc_md = load('/mnt/storage/HaoRan/RRR_motor/M2/ccpc_ll_2.mat');
% ee_md = load('/mnt/storage/rrr_magnum/M2/ccpc_ll_3.mat'); % no baseline but with sd contraint, didn't work too well...
% rsc_md = load('/mnt/storage/HaoRan/RRR_motor/M2/ccpc_ll_3.mat');
% ee_md = load('/mnt/storage/rrr_magnum/M2/ccpc_ll_4.mat'); % with baseline, sd constraint and baseline is lower than amplitude
% rsc_md = load('/mnt/storage/HaoRan/RRR_motor/M2/ccpc_ll_4.mat');

md = cat(2, rsc_md.md2, ee_md.md2);
clusts = cat(2, rsc.clusts{2}, ee.clusts{2});

ratio = [];
% ssize = [];
count = 1;
for ii = 1:length(clusts)
    for jj = 1:length(clusts{ii})
        temp = md{ii}.ll(clusts{ii}{jj}, :);
        ratio{count} = -2 .* (temp(:, 2) - temp(:, 1));
%         ssize{count} = md{ii}.ssize(clusts{ii}{jj});
        count = count + 1;
    end
end

groups = repelem(iscue2, cellfun(@length, ratio));
subratio = cell2mat(ratio');
% ssize = cell2mat(ssize);

figure
hold on
yyaxis left
edges = linspace(min(subratio), max(subratio), 4e2 + 1);
histogram(subratio(~~groups), edges, 'Normalization', 'probability', 'LineStyle', 'none');
ylabel('discrete frequency')
ylim([0 .2])
yyaxis right
cdfplot(subratio(groups == 1));
cdfplot(subratio(groups == 2));
ylabel('cumulative frequency')
xlabel('likelihood ratio (\chi^2)')
% xlim([-600, 1400])
% xlim([-400, 400])
legend('combined', 'cue', 'traj')

subratio = subratio(~~groups);
% ssize = ssize(~~groups);
groups = groups(~~groups);
mu = diff(accumarray(groups, subratio, [2, 1], @mean));

p = zeros(length(groups), 1);
for ii = 1:1e5
    idx = randperm(length(groups));
    p(ii) = diff(accumarray(groups(idx), subratio, [2, 1], @mean));
end

figure
edges = linspace(min(p), max(p), 2e2 + 1);
histogram(p, edges, 'Normalization', 'probability', 'LineStyle', 'none');
xline([0, mu])
% xlim([-40, 40]);
% ylim([0, .02])


%% number of neurons per ensemble
count = cellfun(@length, ratio);
edges = min(count)-.5:max(count)+.5;
counts = zeros(length(edges) - 1, length(unique(iscue2)));
for ii = unique(iscue2(:))'
    counts(:, ii + 1) = histcounts(count(iscue2 == ii), edges);
end

figure
% subplot(1, 2, 1);
bar(edges(1:end-1) + .5, counts(:, 2:3) ./ sum(counts(:, 2:3)), 1);
xlabel('number of neurons');
ylabel('proportion');
ylim([0, .5]);
[p, chi] = chisq2(count(iscue2 == 1), count(iscue2 == 2));
title(['\chi^2 test: p = ' num2str(p) '; \chi^2 = ' num2str(chi)]);

% subplot(1, 2, 2);
% boxplot(count, iscue2);


%% decoding plots
ee_md = load('/mnt/storage/rrr_magnum/M2/ens_bayes.mat');
rsc_md = load('/mnt/storage/HaoRan/RRR_motor/M2/ens_bayes.mat');

md = cat(2, rsc_md.md2, ee_md.md2);
clusts = cat(2, rsc.clusts{2}, ee.clusts{2});

acc = cell2mat(cellfun(@(x) x.acc, md, 'UniformOutput', false));
individual = cell2mat(cellfun(@(x) x.individual, md, 'UniformOutput', false));
err = cell2mat(cellfun(@(x) x.err, md, 'UniformOutput', false));

figure
subplot(1, 2, 1)
temp = err(:, iscue2 == 1);
[~, order] = min(temp);
[~, order] = sort(order);
imagesc(temp(:, order)');
subplot(1, 2, 2)
temp = err(:, iscue2 == 2);
[~, order] = min(temp);
[~, order] = sort(order);
imagesc(temp(:, order)');

figure
boxplot(acc(~~iscue2), iscue2(~~iscue2))
ylim([0 .8])
p = ranksum(acc(iscue2 == 1), acc(iscue2 == 2));
title(['ranksum p = ', num2str(p)])

%% prototyping
clear all

% rest = ensemble('/mnt/md0/Data/RSC_M2/RSC037/2017_09_13/2017_09_13_3.abf');
rest = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_09_13/2017_09_13_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EC002/2018_12_18/2018_12_18_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EC002/2018_12_20/2018_12_20_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EE006/2019_06_05/2019_06_05_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EC005/2019_05_16/2019_05_16_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EC005/2019_05_22/2019_05_22_3.abf');
% rest = ensemble('/mnt/storage/rrr_magnum/M2/EE001/2019_01_07/2019_01_07_3.abf');
rest.set_ops('e_size',5);
rest.set_ops('clust_method','thres');
rest.set_ops('sig', .2);
rest.remove_mvt;
rest.cluster

md = ccpc_ll(rest);
% md = ens_bayes(rest);


%% run batch analysis --- bayes
clear all

traj_thres = .2;

root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));

% root = '/mnt/storage/rrr_magnum/M2';
% animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

md1 = []; md2 = [];
count = 1;
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
%         rest1.hPICA;
        md1{count} = ens_bayes(rest1);
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
%         rest2.hPICA;
        md2{count} = ens_bayes(rest2);
        
        count = count + 1;
    end
end