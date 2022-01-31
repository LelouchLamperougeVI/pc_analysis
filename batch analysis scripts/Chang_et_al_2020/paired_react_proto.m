%% classify them ensembles... again...

clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'pc_list');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'pc_list');

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
hiepi_psth = cat(1, rsc.hiepi_psth{2}, ee.hiepi_psth{2});
session = 1:length(hiepi_psth);
session = repelem(session, cellfun(@length, hiepi_psth));
session = session(:);

hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = cellfun(@(x) mean(x, 2), hiepi_psth, 'UniformOutput', false);
hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = fast_smooth(hiepi_psth, 5 / 150 * 50);
hiepi_psth = (hiepi_psth - min(hiepi_psth, [], 1)) ./ range(hiepi_psth, 1);

hiepi_z = cat(1, rsc.hiepi_z{2}, ee.hiepi_z{2});
hiepi_z = cellfun(@(x) mat2cell(x, size(x, 1), ones(1, size(x, 2))), hiepi_z, 'UniformOutput', false);
hiepi_z = cat(2, hiepi_z{:});
hiepi_z = hiepi_z(:);
% hiepi_z = cellfun(@(x) fast_smooth(x, 4), hiepi_z, 'UniformOutput', false);
hiepi_z = cellfun(@(x) zscore(x(~isnan(x))), hiepi_z, 'UniformOutput', false);


%% conduct paired analysis
aggr = [];
r = [];
for s = unique(session(:))'
    if ~all(ismember(1:2, iscue2(session == s)))
        continue
    end
    
    cue_feat = hiepi_psth(:, iscue2 == 1 & session == s);
    traj_feat = hiepi_psth(:, iscue2 == 2 & session == s);
    rf = corr(cue_feat, traj_feat);
    
    cue_z = hiepi_z(iscue2 == 1 & session == s);
    traj_z = hiepi_z(iscue2 == 2 & session == s);
    rx = nan(length(cue_z), length(traj_z));
    lag = nan(length(cue_z), length(traj_z));
    for ii = 1:length(cue_z)
        for jj = 1:length(traj_z)
            [temp, t] = xcorr(cue_z{ii}, traj_z{jj}, 1 * 19, 'unbiased');
            [rx(ii, jj), idx] = max(temp);
            lag(ii, jj) = t(idx);
            
            r = cat(1, r, temp(:)');
        end
    end
    
    aggr = cat(1, aggr, [lag(:), rx(:), rf(:)]);
end


%% plot xcorr
figure
idx = max(r, [], 2);
[~, idx] = sort(idx);
imagesc(r(idx, :))
colormap jet
xline(ceil(size(r, 2) / 2))


%%
discretize(aggr(:, 3), [-1, -.4, -.2, .2, .4, 1]);

%%
ax(1) = subplot(1, 2, 1);
imagesc(r(idx, :))
colormap jet
ax(2) = subplot(1, 2, 2);
plot(1:size(aggr, 1), aggr(idx, 3))
plot(aggr(idx, 3), 1:size(aggr, 1))
set(ax(1),'YDir','normal')
linkaxes(ax, 'y')



%% hopfield net modelling
thres = .8;
perms = 1e4;

none_ens = find(iscue2 == 0);
cue_ens = find(iscue2 == 1);
traj_ens = find(iscue2 == 2);

test = [randsample(none_ens, perms, true), randsample(traj_ens, perms, true)];

partials = hiepi_psth(:, iscue2 == 1);
partials = double(partials > thres);
partials(~partials) = -1;

% partials = -partials; % interneuron

retrieved = zeros(size(partials, 2), perms);
for ii = 1:perms
    patterns = hiepi_psth(:, test(ii, :));
    patterns = double(patterns > thres);
    patterns(~patterns) = -1;
    
    retrieved(:, ii) = simple_hop(patterns, partials);
end

retrieved(isnan(retrieved)) = 0;
aggr = arrayfun(@(x) histcounts(retrieved(:, x), 0:3) ./ length(retrieved(:, x)), 1:perms, 'UniformOutput', false);
aggr = cell2mat(aggr');

g = repmat({'dnc', 'none', 'traj'}, [perms, 1]);
g = gramm('y', aggr(:), 'x', g);
g.stat_summary('type', '95percentile', 'geom',{'bar','black_errorbar'});
figure
g.draw();

hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue2), 1], @mean);
hi_traj = hi_traj(iscue2 == 2);
hi_traj = traj_ens(hi_traj > prctile(hi_traj, 80));

hi_cue = sum(retrieved == 2, 2);
hi_cue = cue_ens(hi_cue > prctile(hi_cue, 80));

figure
subplot(2, 2, 1)
histogram(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, hi_cue)), 30)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(2, 2, 3)
imagesc(temp(:, idx)')

hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue2), 1], @mean);
hi_traj = hi_traj(iscue2 == 2);
hi_traj = traj_ens(hi_traj < prctile(hi_traj, 20));

hi_cue = sum(retrieved == 2, 2);
hi_cue = cue_ens(hi_cue < prctile(hi_cue, 20));


subplot(2, 2, 2)
histogram(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, hi_cue)), 30)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(2, 2, 4)
imagesc(temp(:, idx)')


%% hopfield paired
aggr = [];
r = [];
retrieved = [];
for s = unique(session(:))'
    if ~all(ismember(1:2, iscue2(session == s)))
        continue
    end
    
    cue_feat = hiepi_psth(:, iscue2 == 1 & session == s);
    traj_feat = hiepi_psth(:, iscue2 == 2 & session == s);
    rf = corr(cue_feat, traj_feat);
    
    cue_z = hiepi_z(iscue2 == 1 & session == s);
    traj_z = hiepi_z(iscue2 == 2 & session == s);
    rx = nan(length(cue_z), length(traj_z));
    lag = nan(length(cue_z), length(traj_z));
    for ii = 1:length(cue_z)
        for jj = 1:length(traj_z)
            [temp, t] = xcorr(cue_z{ii}, traj_z{jj}, 1 * 19, 'unbiased');
            [rx(ii, jj), idx] = max(temp);
            lag(ii, jj) = t(idx);
            
            r = cat(1, r, temp(:)');
        end
    end
    
    aggr = cat(1, aggr, [lag(:), rx(:), rf(:)]);
    
    for ii = 1:size(rx, 1)
        tests = find(rx(ii, :) > .1);
        controls = find(rx(ii, :) < .1);
        if isempty(tests) || isempty(controls)
            continue
        end
        
        partials = cue_feat(:, ii);
        partials = double(partials > thres);
        partials(~partials) = -1;

        for jj = 1:length(tests)
            for kk = 1:length(controls)
                patterns = [traj_feat(:, tests(jj)), traj_feat(:, controls(kk))];
                patterns = double(patterns > thres);
                patterns(~patterns) = -1;

                retrieved = cat(1, retrieved, simple_hop(patterns, partials));
            end
        end
    end
end

p = binopdf(sum(retrieved == 1), sum(~isnan(retrieved)), .5)


