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


%%
animals = regexpi(rsc.file, "RSC...", 'match');
animals = [animals{:}];
temp = regexpi(ee.file, "E[EC]...", 'match');
temp = [temp{:}];
animals = cat(1, animals(:), temp(:));
[a, idx] = unique(animals, 'stable');
idx = cat(1, idx(:), length(animals) + 1);
count = diff(idx);

tab = cat(2, a, num2cell(count));
tab{end + 1, 1} = 'total';
tab{end, 2} = sum(count);

% lfp
bad_lfp = {'RSC036', 'RSC037', 'RSC038', 'EC002', 'EE003'};
col = size(tab, 2) + 1;
for ii = 1:length(a)
    if ismember(a{ii}, bad_lfp)
        tab{ii, col} = 'N/A';
    else
        tab{ii, col} = 'Avail.';
    end
end
tab{end, col} = sum(cellfun(@(x) strcmpi(x, 'avail.'), tab(:, col)));

% cross days
bad_days = {'RSC036', 'RSC037', 'RSC038', 'EC003'};
col = size(tab, 2) + 1;
for ii = 1:length(a)
    if ismember(a{ii}, bad_days)
        tab{ii, col} = 'N/A';
    else
        tab{ii, col} = 'Avail.';
    end
end
tab{end, col} = sum(cellfun(@(x) strcmpi(x, 'avail.'), tab(:, col)));

% rest1 num ensembles: none, cue, traj, all
clusts = cat(1, rsc.clusts{1}(:), ee.clusts{1}(:));
acc = 1; % accumulator
col = size(tab, 2) + 1;
for ii = 1:length(a)
    count = cellfun(@length, clusts(idx(ii):(idx(ii+1) - 1)));
    count = sum(count);
    temp = histcounts(iscue1(acc:(count+acc-1)), (0:3) - .5);
    temp(end+1) = count;
    tab(ii, (1:length(temp)) + col - 1) = num2cell(temp);
    acc = acc + count;
end
tab(end, (1:length(temp)) + col - 1) = num2cell(sum(cell2mat(tab(1:end-1, (1:length(temp)) + col - 1))));

% rest2 num ensembles: none, cue, traj, all
clusts = cat(1, rsc.clusts{2}(:), ee.clusts{2}(:));
acc = 1; % accumulator
col = size(tab, 2) + 1;
for ii = 1:length(a)
    count = cellfun(@length, clusts(idx(ii):(idx(ii+1) - 1)));
    count = sum(count);
    temp = histcounts(iscue2(acc:(count+acc-1)), (0:3) - .5);
    temp(end+1) = count;
    tab(ii, (1:length(temp)) + col - 1) = num2cell(temp);
    acc = acc + count;
end
tab(end, (1:length(temp)) + col - 1) = num2cell(sum(cell2mat(tab(1:end-1, (1:length(temp)) + col - 1))));

cell2table(tab);



