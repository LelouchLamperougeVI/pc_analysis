clear all

root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));

% root = '/mnt/storage/rrr_magnum/M2';
% animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

count = 1;

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    disp(['running ' num2str(a) '/' num2str(length(animals))]);
    
    for s = 1:length(sessions)
        for t = 1:3
            path = fullfile(root, animals{a}, sessions{s}, num2str(t));
            flag = dir(fullfile(path, 'plane1'));
            if isempty(flag)
                path = fullfile(path, 'Plane1');
            else
                path = fullfile(path, 'plane1');
            end
            load(fullfile(path, 'deconv.mat'));
            load(fullfile(path, 'timecourses.mat'));
            
            filt_d = ca_filt(deconv);
            
            dff = tcs.ratio;
            fr{count}(:, t) = sum(deconv > 0) ./ tcs.tt(end);
            fr_d{count}(:, t) = sum(filt_d > 0) ./ tcs.tt(end);
            mu{count}(:, t) = mean(dff);
        end
        disp([num2str(s) '/' num2str(length(sessions))]);
        count = count + 1;
    end
end

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    disp(['running ' num2str(a) '/' num2str(length(animals))]);
    
    for s = 1:length(sessions)
        for t = 1:3
            path = fullfile(root, animals{a}, sessions{s}, num2str(t));
            flag = dir(fullfile(path, 'plane1'));
            if isempty(flag)
                path = fullfile(path, 'Plane1');
            else
                path = fullfile(path, 'plane1');
            end
            load(fullfile(path, 'deconv.mat'));
            load(fullfile(path, 'timecourses.mat'));
            
            filt_d = ca_filt(deconv);
            
            dff = tcs.ratio;
            fr{count}(:, t) = sum(deconv > 0) ./ tcs.tt(end);
            fr_d{count}(:, t) = sum(filt_d > 0) ./ tcs.tt(end);
            mu{count}(:, t) = mean(dff);
        end
        disp([num2str(s) '/' num2str(length(sessions))]);
        count = count + 1;
    end
end


% arrayfun(@(ii) corr(mu{ii}(:, 2), fr_d{ii}(:, 2)), 1:length(mu));


%%
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack', 'hiepi_struct');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack', 'hiepi_struct');

% classify ensembles
clust_stacks{1} = cat(1, rsc.clust_stacks{1}, ee.clust_stacks{1});
clust_stacks{2} = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});
% clust_stacks{1} = rsc.clust_stacks{1};
% clust_stacks{2} = rsc.clust_stacks{2};
% clust_stacks{1} = ee.clust_stacks{1};
% clust_stacks{2} = ee.clust_stacks{2};

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

clusts1 = cat(2, rsc.clusts{1}, ee.clusts{1});
clusts2 = cat(2, rsc.clusts{2}, ee.clusts{2});
hiepi_struct1 = cat(1, rsc.hiepi_struct{1}{:}, ee.hiepi_struct{1}{:});
hiepi_struct2 = cat(1, rsc.hiepi_struct{2}{:}, ee.hiepi_struct{2}{:});

clearvars -except rsc ee iscue1 iscue2 clusts1 clusts2 hiepi_struct1 hiepi_struct2

%% REST metrics
load('review_extra_metrics.mat')

mu_fr = [];
mu_dff = [];
for ii = 1:length(clusts1)
    for jj = 1:length(clusts1{ii})
        mu_fr = [mu_fr; fr_d{ii}(clusts1{ii}{jj}, 1)];
        mu_dff = [mu_dff; mu{ii}(clusts1{ii}{jj}, 1)];
    end
end

for ii = 1:length(clusts2)
    for jj = 1:length(clusts2{ii})
        mu_fr = [mu_fr; fr_d{ii}(clusts2{ii}{jj}, 3)];
        mu_dff = [mu_dff; mu{ii}(clusts2{ii}{jj}, 3)];
    end
end

mu_dff = double(mu_dff);

idx = cat(2, cat(2, clusts1{:}), cat(2, clusts2{:}));
iscue = repelem([iscue1; iscue2], cellfun(@length, idx));

figure
subplot(2, 2, 1)
cdfplot(mu_dff(iscue == 1))
hold on
cdfplot(mu_dff(iscue == 2))
[~, p] = kstest2(mu_dff(iscue == 1), mu_dff(iscue == 2));
xlim([0 30])
xlabel('mean \Delta F/F (%)')
title(['kstest2 p = ' num2str(p)])

subplot(2, 2, 2)
boxplot(mu_dff(iscue ~= 0), iscue(iscue ~= 0))
p = ranksum(mu_dff(iscue == 1), mu_dff(iscue == 2));
ylim([0 20])
ylabel('mean \Delta F/F (%)')
title(['ranksum p = ' num2str(p)])

subplot(2, 2, 3)
cdfplot(mu_fr(iscue == 1))
hold on
cdfplot(mu_fr(iscue == 2))
[~, p] = kstest2(mu_fr(iscue == 1), mu_fr(iscue == 2));
xlim([0 2])
xlabel('mean calcium spike rate (events/sec)')
title(['kstest2 p = ' num2str(p)])

subplot(2, 2, 4)
boxplot(mu_fr(iscue ~= 0), iscue(iscue ~= 0))
p = ranksum(mu_fr(iscue == 1), mu_fr(iscue == 2));
ylim([.2 1.4])
ylabel('mean calcium spike rate (events/sec)')
title(['ranksum p = ' num2str(p)])


%% RUN metrics
mu_fr = [];
mu_dff = [];
for ii = 1:length(clusts1)
    for jj = 1:length(clusts1{ii})
        mu_fr = [mu_fr; fr_d{ii}(clusts1{ii}{jj}, 2)];
        mu_dff = [mu_dff; mu{ii}(clusts1{ii}{jj}, 2)];
    end
end

for ii = 1:length(clusts2)
    for jj = 1:length(clusts2{ii})
        mu_fr = [mu_fr; fr_d{ii}(clusts2{ii}{jj}, 2)];
        mu_dff = [mu_dff; mu{ii}(clusts2{ii}{jj}, 2)];
    end
end

mu_dff = double(mu_dff);

idx = cat(2, cat(2, clusts1{:}), cat(2, clusts2{:}));
iscue = repelem([iscue1; iscue2], cellfun(@length, idx));

figure
subplot(2, 2, 1)
cdfplot(mu_dff(iscue == 1))
hold on
cdfplot(mu_dff(iscue == 2))
[~, p] = kstest2(mu_dff(iscue == 1), mu_dff(iscue == 2));
xlim([0 60])
xlabel('mean \Delta F/F (%)')
title(['kstest2 p = ' num2str(p)])

subplot(2, 2, 2)
boxplot(mu_dff(iscue ~= 0), iscue(iscue ~= 0))
p = ranksum(mu_dff(iscue == 1), mu_dff(iscue == 2));
ylim([0 40])
ylabel('mean \Delta F/F (%)')
title(['ranksum p = ' num2str(p)])

subplot(2, 2, 3)
cdfplot(mu_fr(iscue == 1))
hold on
cdfplot(mu_fr(iscue == 2))
[~, p] = kstest2(mu_fr(iscue == 1), mu_fr(iscue == 2));
xlim([0 2])
xlabel('mean calcium spike rate (events/sec)')
title(['kstest2 p = ' num2str(p)])

subplot(2, 2, 4)
boxplot(mu_fr(iscue ~= 0), iscue(iscue ~= 0))
p = ranksum(mu_fr(iscue == 1), mu_fr(iscue == 2));
ylim([.2 1.2])
ylabel('mean calcium spike rate (events/sec)')
title(['ranksum p = ' num2str(p)])


%%
iscue = cat(1, iscue1, iscue2);
rest = repelem((1:2)', [length(iscue1); length(iscue2)]);
hiepi_struct1(56).on = nan;
ISI = cat(1, arrayfun(@(x) mean(diff(x.on)), hiepi_struct1), arrayfun(@(x) mean(diff(x.on)), hiepi_struct2));
rate = cat(1, arrayfun(@(x) length(x.on) / (x.on(length(x.on)) - x.on(1)), hiepi_struct1), arrayfun(@(x) length(x.on) / (x.on(length(x.on)) - x.on(1)), hiepi_struct2));

figure
subplot(1, 4, 1)
boxplot(ISI(iscue ~= 0 & rest == 1), iscue(iscue ~= 0 & rest == 1))
p = ranksum(ISI(iscue == 1 & rest == 1), ISI(iscue == 2 & rest == 1));
ylim([0 50])
ylabel('inter-event-interval (sec)')
title(['ranksum p = ' num2str(p)])

subplot(1, 4, 2)
boxplot(ISI(iscue ~= 0 & rest == 2), iscue(iscue ~= 0 & rest == 2))
p = ranksum(ISI(iscue == 1 & rest == 2), ISI(iscue == 2 & rest == 2));
ylim([0 50])
ylabel('inter-event-interval (sec)')
title(['ranksum p = ' num2str(p)])

subplot(1, 4, 3)
boxplot(rate(iscue ~= 0 & rest == 1), iscue(iscue ~= 0 & rest == 1))
p = ranksum(rate(iscue == 1 & rest == 1), rate(iscue == 2 & rest == 1));
ylim([0 .2])
ylabel('reactivation rate (events/sec)')
title(['ranksum p = ' num2str(p)])

subplot(1, 4, 4)
boxplot(rate(iscue ~= 0 & rest == 2), iscue(iscue ~= 0 & rest == 2))
p = ranksum(rate(iscue == 1 & rest == 2), rate(iscue == 2 & rest == 2));
ylim([0 .2])
ylabel('reactivation rate (events/sec)')
title(['ranksum p = ' num2str(p)])