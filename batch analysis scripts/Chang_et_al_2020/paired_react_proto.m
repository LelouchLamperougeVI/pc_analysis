%% classify them ensembles... again...

clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'hiepi_struct', 'pc_list', 'file');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'hiepi_z', 'hiepi_psth', 'hiepi_lfp_pw', 'hiepi_struct', 'pc_list', 'file');

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




%% conduct paired analysis
if rest == 1
    iscue = iscue1;
elseif rest == 2
    iscue = iscue2;
end

aggr = [];
r = [];
for s = unique(session(:))'
    if ~all(ismember(1:2, iscue(session == s)))
        continue
    end
    
    cue_feat = hiepi_psth(:, iscue == 1 & session == s);
    traj_feat = hiepi_psth(:, iscue == 2 & session == s);
    rf = corr(cue_feat, traj_feat);
    
    cue_z = hiepi_z(iscue == 1 & session == s);
    traj_z = hiepi_z(iscue == 2 & session == s);
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

figure
subplot(1, 2, 1)
idx = max(r, [], 2);
[~, idx] = sort(idx);
imagesc(r(idx, :))
colormap jet
xline(ceil(size(r, 2) / 2))

subplot(1, 2, 2)
cdfplot(aggr(:, 2))
knee = prctile(aggr(:, 2), 80);
xline(knee);
title(['elbow of curve at ~80 percentile = ' num2str(knee)])



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
thres = .5;
perms = 1e4;

none_ens = find(iscue == 0);
cue_ens = find(iscue == 1);
traj_ens = find(iscue == 2);

test = [randsample(traj_ens, perms, true), randsample(traj_ens, perms, true)];

partials = hiepi_psth(:, iscue == 1);
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

figure
histogram(diff(aggr(:, 2:3), 1, 2), 50)
xline(prctile(diff(aggr(:, 2:3), 1, 2), [25, 50, 75]));

figure
hi_traj = traj_ens;
[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]);
hi_cue = hi_cue(:, hi_traj);
[~, hi_cue] = max(hi_cue);
subplot(1, 3, 1)
histogram(diag(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, cue_ens(hi_cue)))), 10)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(1, 3, 2)
imagesc(temp(:, idx)')
temp = hiepi_psth(:, cue_ens(hi_cue));
subplot(1, 3, 3)
imagesc(temp(:, idx)')

figure
hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue), 1], @mean);
hi_traj = hi_traj(iscue == 2);
hi_traj = traj_ens(hi_traj > prctile(hi_traj, 75));
[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]);
hi_cue = hi_cue(:, hi_traj);
[~, hi_cue] = max(hi_cue);
hi_cue = 1:length(cue_ens);
subplot(3, 4, 1)
histogram(diag(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, cue_ens(hi_cue)))), 10)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 5)
imagesc(temp(:, idx)')
temp = hiepi_psth(:, cue_ens(hi_cue));
subplot(3, 4, 9)
imagesc(temp(:, idx)')

hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue), 1], @mean);
hi_traj = hi_traj(iscue == 2);
hi_traj = traj_ens(hi_traj > prctile(hi_traj, 50) & hi_traj <= prctile(hi_traj, 75));
[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]);
hi_cue = hi_cue(:, hi_traj);
[~, hi_cue] = max(hi_cue);
hi_cue = 1:length(cue_ens);
subplot(3, 4, 2)
histogram(diag(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, cue_ens(hi_cue)))), 10)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 6)
imagesc(temp(:, idx)')
temp = hiepi_psth(:, cue_ens(hi_cue));
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 10)
imagesc(temp(:, idx)')

hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue), 1], @mean);
hi_traj = hi_traj(iscue == 2);
hi_traj = traj_ens(hi_traj > prctile(hi_traj, 25) & hi_traj <= prctile(hi_traj, 50));
[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]);
hi_cue = hi_cue(:, hi_traj);
[~, hi_cue] = max(hi_cue);
hi_cue = 1:length(cue_ens);
subplot(3, 4, 3)
histogram(diag(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, cue_ens(hi_cue)))), 10)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 7)
imagesc(temp(:, idx)')
temp = hiepi_psth(:, cue_ens(hi_cue));
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 11)
imagesc(temp(:, idx)')

hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue), 1], @mean);
hi_traj = hi_traj(iscue == 2);
hi_traj = traj_ens(hi_traj <= prctile(hi_traj, 25));
[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]);
hi_cue = hi_cue(:, hi_traj);
[~, hi_cue] = max(hi_cue);
hi_cue = 1:length(cue_ens);
subplot(3, 4, 4)
histogram(diag(corr(hiepi_psth(:, hi_traj), hiepi_psth(:, cue_ens(hi_cue)))), 10)
temp = hiepi_psth(:, hi_traj);
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 8)
imagesc(temp(:, idx)')
temp = hiepi_psth(:, cue_ens(hi_cue));
[~, idx] = max(temp);
[~, idx] = sort(idx);
subplot(3, 4, 12)
imagesc(temp(:, idx)')


%%
hi_traj = traj_ens;
temp = hiepi_psth(:, hi_traj);
ham = permute(temp > thres, [2 1]) ~= permute(hiepi_psth(:, cue_ens) > thres, [3 1 2]);
ham = squeeze(sum(ham, 2))' ./ size(hiepi_psth, 1);

[X, Y] = find(retrieved == 2);
hi_cue = accumarray([X, test(Y, 2)], 1, [length(cue_ens), size(hiepi_psth, 2)]) ./ accumarray(test(:, 2), 1, [size(hiepi_psth, 2), 1])';
success = hi_cue(:, traj_ens);

g = gramm('x', ham, 'y', success);
g.geom_jitter('width', .02, 'height', .01);
g.set_names('x', 'normalized hamming distance', 'y', 'fraction of successful retrievals');
figure
g.draw;

[X, Y] = find(ham < .3);
temp = hiepi_psth(:, traj_ens(Y));
figure
subplot(1, 2, 1);
[~, idx] = max(temp);
[~, idx] = sort(idx);
imagesc(temp(:, idx)');
subplot(1, 2, 2);
temp = hiepi_psth(:, cue_ens(X));
imagesc(temp(:, idx)');

[X, Y] = find(success > .75);
ol = @(x) sum(hiepi_psth(:, cue_ens(X(x))) > thres & hiepi_psth(:, traj_ens(Y(x))) > thres, 2) ./ sum(hiepi_psth(:, traj_ens(Y(x))) > thres, 2);
overlap = ol(1:length(Y));
temp = arrayfun(@(x) ol(datasample(1:length(Y), length(Y))), 1:1e3, 'UniformOutput', false);
ci = cell2mat(temp);
% ci = prctile(ci, [2.5 97.5], 2);

[X, Y] = find(success < .25);
ol = @(x) sum(hiepi_psth(:, cue_ens(X(x))) > thres & hiepi_psth(:, traj_ens(Y(x))) > thres, 2) ./ sum(hiepi_psth(:, traj_ens(Y(x))) > thres, 2);
overlap = cat(1, overlap, ol(1:length(Y)));
temp = arrayfun(@(x) ol(datasample(1:length(Y), length(Y))), 1:1e3, 'UniformOutput', false);
ci = cat(1, ci, cell2mat(temp));
ci = prctile(ci, [2.5 97.5], 2);

x = repmat((1:50)', 2, 2);
y = cat(2, overlap(:), 1 - overlap(:));
facet = repmat({'overlap', 'no overlap'}, length(overlap), 1);
ymin = cat(2, ci(:, 1), 1 - ci(:, 2));
ymax = cat(2, ci(:, 2), 1 - ci(:, 1));
color = repmat(repelem({'>75%'; '<25%'}, [length(overlap) / 2, length(overlap) / 2]), 1, 2);
g = gramm('x', x(:), 'y', y(:), 'color', color(:), 'ymin', ymin(:), 'ymax', ymax(:));
g.facet_grid([], facet(:));
g.geom_line();
g.geom_interval();
g.axe_property('XLim', [1 50], 'YLim', [0 1]);
g.set_names('x', 'position', 'y', 'fraction overlapping pixels', 'color', 'retrieval success rate');
figure
g.draw();

[X, Y] = find(success > .75);
temp = corr(hiepi_psth(:, cue_ens(X)), hiepi_psth(:, traj_ens(Y)));
r = diag(temp);
g = repelem(1, length(Y));
[X, Y] = find(success < .25);
temp = corr(hiepi_psth(:, cue_ens(X)), hiepi_psth(:, traj_ens(Y)));
r = cat(1, r, diag(temp));
g = cat(2, g, repelem(2, length(Y)));
figure
boxplot(r, g);
p = ranksum(r(g == 1), r(g == 2));
title(['ranksum p-value = ' num2str(p)])

[X, Y] = find(success > .75);
ol = @(x) sum(hiepi_psth(:, cue_ens(X(x))) > thres & hiepi_psth(:, traj_ens(Y(x))) > thres, 1) ./ sum(hiepi_psth(:, traj_ens(Y(x))) > thres, 1);
overlap = ol(1:length(Y));
g = repelem(1, length(Y));
[X, Y] = find(success < .25);
ol = @(x) sum(hiepi_psth(:, cue_ens(X(x))) > thres & hiepi_psth(:, traj_ens(Y(x))) > thres, 1) ./ sum(hiepi_psth(:, traj_ens(Y(x))) > thres, 1);
overlap = cat(2, overlap, ol(1:length(Y)));
g = cat(2, g, repelem(2, length(Y)));
figure
boxplot(overlap, g);
p = ranksum(overlap(g == 1), overlap(g == 2));
title(['ranksum p-value = ' num2str(p)])


%% hopfield paired
knee = 0.061289393424153;
aggr = [];
r = [];
retrieved = [];
ham = [];
for s = unique(session(:))'
    if ~all(ismember(1:2, iscue(session == s)))
        continue
    end
    
    cue_feat = hiepi_psth(:, iscue == 1 & session == s);
    traj_feat = hiepi_psth(:, iscue == 2 & session == s);
    rf = corr(cue_feat, traj_feat);
    
    cue_idx = find(iscue == 1 & session == s);
    traj_idx = find(iscue == 2 & session == s);
    cue_z = hiepi_z(cue_idx);
    traj_z = hiepi_z(traj_idx);
    [traj_idx, cue_idx] = meshgrid(traj_idx, cue_idx);
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
    
    aggr = cat(1, aggr, [lag(:), rx(:), rf(:), cue_idx(:), traj_idx(:)]);
    
    for ii = 1:size(rx, 1)
        tests = find(rx(ii, :) > knee);
        controls = setxor(1:length(rx(ii, :)), tests);
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
                ham = cat(1, ham, sum(patterns ~= partials));
            end
        end
    end
end
ham = ham ./ size(hiepi_psth, 1);

p = binopdf(sum(retrieved == 1), sum(~isnan(retrieved)), .5);
figure
bar(histcounts(retrieved) ./ sum(~isnan(retrieved)))
title(['binomial test p-value = ' num2str(p)])

idx = max(r, [], 'all'); % remove outlier
[idx, ~] = find(r == idx);
r(idx, :) = [];

figure
subplot(2, 2, 3)
plot(mean(r))
xline(ceil(size(r, 2) / 2))
yline(0)
auc = cumtrapz(mean(r)) ./ trapz(mean(r));
title(['AUC = ' num2str(auc(ceil(size(r, 2) / 2)))])

subplot(2, 2, 2)
cdfplot(aggr(:, 2))
knee = prctile(aggr(:, 2), 80);
xline(knee);
title(['elbow of curve at ~80 percentile = ' num2str(knee)])

subplot(2, 2, 4)
boxplot(aggr(aggr(:, 2) > knee, 1))
p = signrank(aggr(aggr(:, 2) > knee, 1));
title(['signrank p-value = ' num2str(p)])
ylim([-19 19])

subplot(2, 2, 1)
idx = max(r, [], 2);
[temp, idx] = sort(idx);
imagesc(r(idx, :))
colormap viridis
xline(ceil(size(r, 2) / 2))
yline(find(temp > knee, 1))

figure
subplot(2, 2, 1)
temp = hiepi_psth(:, aggr(aggr(:, 2) > knee, 4));
[~, idx] = max(temp);
[~, idx] = sort(idx);
imagesc(temp(:, idx)')
subplot(2, 2, 2)
temp = hiepi_psth(:, aggr(aggr(:, 2) > knee, 5));
imagesc(temp(:, idx)')

subplot(2, 2, 3)
temp = hiepi_psth(:, aggr(aggr(:, 2) <= knee, 4));
[~, idx] = max(temp);
[~, idx] = sort(idx);
imagesc(temp(:, idx)')
subplot(2, 2, 4)
temp = hiepi_psth(:, aggr(aggr(:, 2) <= knee, 5));
imagesc(temp(:, idx)')

figure
g = (aggr(:, 2) > knee) + 1;
boxplot(aggr(:, 3), g);
p = ranksum(aggr(aggr(:, 2) > knee, 3), aggr(aggr(:, 2) <= knee, 3));
title(['ranksum p-value = ' num2str(p)])

figure
boxplot(ham)
p = signrank(diff(ham, 1, 2));
title(['signrank p-value = ' num2str(p)])
ylim([0 1])


%% xcorr shuffle instead of knee
shuffles = 1e3;
r = [];
r_shuffle = [];
for ii = 1:size(aggr, 1)
    cue_z = zscore(hiepi_z{aggr(ii, 4)});
    traj_z = zscore(hiepi_z{aggr(ii, 5)});
    temp = xcorr(cue_z, traj_z, 1 * 19, 'unbiased');
    r = cat(1, r, temp(:)');
    
    k = randsample(length(cue_z), shuffles, true);
    s_stack = zeros(shuffles, length(temp));
    parfor s = 1:shuffles
        temp = circshift(cue_z, k(s));
        s_stack(s, :) = xcorr(temp, traj_z, 1 * 19, 'unbiased');
    end
    temp = [mean(s_stack); prctile(s_stack, 99); sum(r(end, :) > s_stack, 1) ./ shuffles];
    r_shuffle = cat(1, r_shuffle, permute(temp, [3 2 1]));
end

couple_idx = max(r, [], 2) > max(r_shuffle(:, :, 2), [], 2);

% t = linspace(-1, 1, size(r, 2));
% t = repmat(t, size(r, 1) + size(r_shuffle, 1), 1);
% color = repelem({'real'; 'shuffle'}, [size(r, 1); size(r_shuffle, 1)], 1);
% g = gramm('x', t, 'y', cat(1, r, r_shuffle(:, :, 1)), 'color', color);
% g.stat_summary('type', 'sem', 'geom', 'area');
% figure
% g.draw();
% 
% p = 1 - r_shuffle(:, :, 3);
% p = p < .05;
% heads = get_head(p');
% tails = get_head(p(:, end:-1:1)'); tails = tails(end:-1:1, :);
% [heads, idx] = find(heads);
% [tails, ~] = find(tails);
% idx = accumarray(idx, tails - heads, [size(p, 1), 1], @max);
% couple_idx = idx >= 10;


%% Detection of react. event onsets and coupling
% edges = linspace(-5, 5, 101 + 1);
% r = nan(size(aggr, 1), 101);
r = [];
for ii = 1:size(aggr, 1)
%     temp = hiepi_struct(aggr(ii, 4)).peaks - hiepi_struct(aggr(ii, 5)).peaks';
%     [~, lags, temp] = point_xcorr(hiepi_struct(aggr(ii, 4)).on, hiepi_struct(aggr(ii, 5)).on, length(hiepi_z{aggr(ii, 4)}), 19);
    [~, lags, temp] = point_xcorr2(hiepi_struct(aggr(ii, 4)).on, hiepi_struct(aggr(ii, 5)).on, isnan(hiepi_z{aggr(ii, 4)}), ts{65}, .1);
    if numel(temp) < 100
        continue
    end
%     r(ii, :) = histcounts(temp, edges) ./ numel(temp);
    idx = find(lags == 0);
    r = cat(1, r, temp(idx-19:idx+19));
end

thres = .15;
frac2 = cellfun(@(x) sum(abs(x) < thres) / numel(x), delay);
dur2 = arrayfun(@(x) dur2{x}(abs(delay{x}) < thres), 1:length(dur2), 'uniformoutput', false);
dur2 = cellfun(@mean, dur2);
delay = cellfun(@(x) mean(x(abs(x) < thres), 'omitnan'), delay);



%%
ol = @(x) sum(hiepi_psth(:, aggr(x, 4)) > thres & hiepi_psth(:, aggr(x, 5)) > thres, 2) ...
    ./ sum(hiepi_psth(:, aggr(x, 5)) > thres, 2);
overlap = ol(aggr(:, 2) > knee);
temp = arrayfun(@(x) ol(datasample(find(aggr(:, 2) > knee), sum(aggr(:, 2) > knee))), 1:1e3, 'UniformOutput', false);
ci = cell2mat(temp);
% ci = prctile(ci, [2.5 97.5], 2);

overlap = cat(1, overlap, ol(aggr(:, 2) <= knee));
temp = arrayfun(@(x) ol(datasample(find(aggr(:, 2) <= knee), sum(aggr(:, 2) <= knee))), 1:1e3, 'UniformOutput', false);
ci = cat(1, ci, cell2mat(temp));
ci = prctile(ci, [2.5 97.5], 2);

x = repmat((1:50)', 2, 2);
y = cat(2, overlap(:), 1 - overlap(:));
facet = repmat({'overlap', 'no overlap'}, length(overlap), 1);
ymin = cat(2, ci(:, 1), 1 - ci(:, 2));
ymax = cat(2, ci(:, 2), 1 - ci(:, 1));
color = repmat(repelem({'>75%'; '<25%'}, [length(overlap) / 2, length(overlap) / 2]), 1, 2);
g = gramm('x', x(:), 'y', y(:), 'color', color(:), 'ymin', ymin(:), 'ymax', ymax(:));
g.facet_grid([], facet(:));
g.geom_line();
g.geom_interval();
g.axe_property('XLim', [1 50], 'YLim', [0 1]);
g.set_names('x', 'position', 'y', 'fraction overlapping pixels', 'color', 'retrieval success rate');
figure
g.draw();


ol = @(x) sum(hiepi_psth(:, aggr(x, 4)) > thres & hiepi_psth(:, aggr(x, 5)) > thres, 1) ...
    ./ sum(hiepi_psth(:, aggr(x, 5)) > thres, 1);
overlap = ol(aggr(:, 2) > knee);
g = repelem(1, sum(aggr(:, 2) > knee));

overlap = cat(2, overlap, ol(aggr(:, 2) <= knee));
g = cat(2, g, repelem(2, sum(aggr(:, 2) <= knee)));
figure
boxplot(overlap, g);
p = ranksum(overlap(g == 1), overlap(g == 2));
title(['ranksum p-value = ' num2str(p)])

%%
r_coup = corr(hiepi_psth(:, aggr(aggr(:, 2) > knee, 4))', hiepi_psth(:, aggr(aggr(:, 2) > knee, 5))');
r_unco = corr(hiepi_psth(:, aggr(aggr(:, 2) <= knee, 4))', hiepi_psth(:, aggr(aggr(:, 2) <= knee, 5))');

figure
subplot(1, 2, 1)
imagesc(r_coup);
axis square
rbmap('caxis', [-1 1], 'interp', 255);
subplot(1, 2, 2)
imagesc(r_unco);
axis square
rbmap('caxis', [-1 1], 'interp', 255);


%%
idx = find(aggr(:, 2) > knee);
idx = cat(2, {idx}, arrayfun(@(x) datasample(idx, length(idx)), 1:1e3, 'UniformOutput', false));
idx = cell2mat(idx);

r_coup = zeros(size(hiepi_psth, 1), size(idx, 2));
for s = 1:size(idx, 2)
    temp = corr(hiepi_psth(:, aggr(idx(:, s), 4))', hiepi_psth(:, aggr(idx(:, s), 5))');
    r_coup(:, s) = diag(temp);
end
ci = prctile(r_coup, [2.5 97.5], 2);

g = gramm('x', 1:50, 'y', r_coup(:, 1), 'ymin', ci(:, 1), 'ymax', ci(:, 2));
g.geom_line();
g.geom_interval();
% g.axe_property('XLim', [1 50], 'YLim', [0 1]);
g.set_names('x', 'position', 'y', 'fraction overlapping pixels', 'color', 'retrieval success rate');
figure
g.draw();

%% mutual information
bins = 10;

idx = find(aggr(:, 2) > knee);
idx = cat(2, {idx}, arrayfun(@(x) datasample(idx, length(idx)), 1:1e3, 'UniformOutput', false));
idx = cell2mat(idx);

mi = zeros(size(hiepi_psth, 1), size(idx, 2));
for s = 1:size(idx, 2)
    X = hiepi_psth(:, aggr(idx(:, s), 4));
    Y = hiepi_psth(:, aggr(idx(:, s), 5));
    X = discretize(X, linspace(0, 1, bins + 1));
    Y = discretize(Y, linspace(0, 1, bins + 1));
    XY = zeros(size(X, 1), bins, bins);
    for ii = 1:size(X, 1)
        XY(ii, :, :) = accumarray([X(ii, :)', Y(ii, :)'], 1, [bins, bins]);
    end
    XY = XY ./ sum(XY, [2, 3]);

    X = arrayfun(@(x) histcounts(X(x, :), 1:(bins+1), 'Normalization', 'probability'), 1:size(X, 1), 'UniformOutput', false);
    X = cell2mat(X');
    Y = arrayfun(@(x) histcounts(Y(x, :), 1:(bins+1), 'Normalization', 'probability'), 1:size(Y, 1), 'UniformOutput', false);
    Y = cell2mat(Y');
    Y = permute(Y, [1 3 2]);

    temp = XY .* log2(XY ./ (X .* Y));
    temp = sum(temp, [2, 3], 'omitnan');
    
    mi(:, s) = temp;
end

ci = prctile(mi, [2.5 97.5], 2);

g = gramm('x', 1:50, 'y', mi(:, 1), 'ymin', ci(:, 1), 'ymax', ci(:, 2));
g.geom_line();
g.geom_interval();
% g.axe_property('XLim', [1 50], 'YLim', [0 1]);
g.set_names('x', 'position', 'y', 'fraction overlapping pixels', 'color', 'retrieval success rate');
figure
g.draw();


X = hiepi_psth(:, aggr(aggr(:, 2) <= knee, 4));
Y = hiepi_psth(:, aggr(aggr(:, 2) <= knee, 5));
X = discretize(X, linspace(0, 1, bins + 1));
Y = discretize(Y, linspace(0, 1, bins + 1));
XY = zeros(size(X, 1), bins, bins);
for ii = 1:size(X, 1)
    XY(ii, :, :) = accumarray([X(ii, :)', Y(ii, :)'], 1, [bins, bins]);
end
XY = XY ./ sum(XY, [2, 3]);

X = arrayfun(@(x) histcounts(X(x, :), 1:(bins+1), 'Normalization', 'probability'), 1:size(X, 1), 'UniformOutput', false);
X = cell2mat(X');
Y = arrayfun(@(x) histcounts(Y(x, :), 1:(bins+1), 'Normalization', 'probability'), 1:size(Y, 1), 'UniformOutput', false);
Y = cell2mat(Y');
Y = permute(Y, [1 3 2]);

mi = XY .* log2(XY ./ (X .* Y));
mi = sum(mi, [2, 3], 'omitnan');

plot(mi)


%% plot pairs
targets = [aggr(aggr(:, 2) > knee, 4), aggr(aggr(:, 2) > knee, 5), aggr(aggr(:, 2) > knee, 3)];

for ii = 1:size(targets, 1)
    if ~mod(ii - 1, 25)
        figure
    end
    
    ax(1) = subplot(15, 5, mod(ii - 1, 25) + 1 + floor(mod(ii-1, 25) / 5) * 10);
    temp = clust_stacks{rest}{targets(ii, 1)}; [~, idx] = max(temp); [~, idx] = sort(idx);
    imagesc(-temp(:, idx)');
    colormap bone
    
    ax(2) = subplot(15, 5, mod(ii - 1, 25) + 1 + floor(mod(ii-1, 25) / 5) * 10 + 5);
    temp = clust_stacks{rest}{targets(ii, 2)}; [~, idx] = max(temp); [~, idx] = sort(idx);
    imagesc(-temp(:, idx)');
    colormap bone
    
    ax(3) = subplot(15, 5, mod(ii - 1, 25) + 1 + floor(mod(ii-1, 25) / 5) * 10 + 10);
    plot(hiepi_psth(:, targets(ii, 1)))
    hold on
    plot(hiepi_psth(:, targets(ii, 2)))
    title(num2str(targets(ii, 3)))
    
    linkaxes(ax, 'x');
    xlim([1 50])
end


%% cue/traj polar plots
temp = cat(1, rsc.hiepi_psth{rest}, ee.hiepi_psth{rest});
temp = cat(2, temp{:});
temp = cellfun(@(x) mean(x, 2), temp, 'UniformOutput', false);
temp = cat(2, temp{:});
temp = repmat(temp, 3, 1); % circular-ish smoothing
temp = fast_smooth(temp, 5 / 150 * 50);
temp = temp((size(temp, 1)/3 + 1):(size(temp, 1)/3*2), :);
temp = (temp - min(temp, [], 1)) ./ range(temp, 1);

[~, match] = arrayfun(@(x) find(belt_num == x), unique(belt_num(belt_num > 0)), 'UniformOutput', false);
% idx = aggr(:, 2) > knee;
% idx = aggr(:, 2) <= knee;
idx = couple_idx;
% idx = ~couple_idx;
targets = hiepi_psth(:, aggr(idx, 4));
[~, targets] = max(targets);

d = @(x, y) min(abs(cat(3, x(:) + size(hiepi_psth, 1) - y(:)', x(:) - y(:)')), [], [2, 3]);

targets = cellfun(@(x) d(targets, x), match, 'UniformOutput', false);
targets = cell2mat(targets);
[~, targets] = min(targets, [], 2);

targets(targets == 1) = 2;
targets = targets - 1;
bn = belt_num;
bn(bn == 1) = 2;
bn = bn - 1;

figure
for ii = 1:length(unique(targets))
    subplot(2, length(unique(targets)), ii);
    
    theta = linspace(0, 2*pi, size(temp, 1) + 1);
    
    polarhistogram(theta(find(bn == ii)), theta);
    hold on
    
    rho = temp(:, aggr(idx, 4));
    rho = mean(rho(:, targets == ii), 2);
    rho = rho ./ max(rho);
    polarplot(theta, [rho; rho(1)]);
    
    rho = temp(:, aggr(idx, 5));
    p = arrayfun(@(x) circ_mean(theta(1:(end-1))', rho(:, x)), find(targets == ii));
    rho = mean(rho(:, targets == ii), 2);
    rho = rho ./ max(rho);
    polarplot(theta, [rho; rho(1)]);
    
    r = circ_r(p(:));
    [phi, ul, ll] = circ_mean(p(:));
    polarplot([0, phi], [0, r]);
    polarplot([0, ul], [0, r], 'k--');
    polarplot([0, ll], [0, r], 'k--');
    
    pval = circ_rtest(p);
    title(['p-value = ' num2str(pval)]);
    
    
    subplot(2, length(unique(targets)), length(unique(targets)) + ii);
    polarhistogram(p, 36);
    title(['n = ' num2str(length(p))])
end

figure
r = diag(corr(temp(:, aggr(idx, 4)), temp(:, aggr(idx, 5))));
boxplot(r, targets)
p = arrayfun(@(x) signrank(r(targets == x)), unique(targets)) .* length(unique(targets));
title(['p-values: ' num2str(p(:)')])


%%
bins = 50;

% hi_traj = accumarray(test(:, 2), aggr(:, 3), [length(iscue2), 1], @mean);
% hi_traj = hi_traj(iscue2 == 2);
% hi_traj = traj_ens(hi_traj > prctile(hi_traj, 75));
% hi_traj = traj_ens(hi_traj < prctile(hi_traj, 25));

% hi_traj = traj_ens(mean(success) > prctile(mean(success), 75));
% hi_traj = traj_ens(mean(success) < prctile(mean(success), 25));

edges = linspace(0, 150, bins+1);
P_se = accumarray([discretize(cell2mat(s2(aggr(aggr(:, 2) > knee, 5))), edges) discretize(cell2mat(e2(aggr(aggr(:, 2) > knee, 5))), edges)], 1, [bins bins]);
% P_se = accumarray([discretize(cell2mat(s2(hi_traj)), edges) discretize(cell2mat(e2(hi_traj)), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_e_cond_s = P_se ./ sum(P_se, 2);

P_e_cond_s(isnan(P_e_cond_s)) = 0;
P_e_cond_s = imgaussfilt(P_e_cond_s, 2);

figure
subplot(2, 2, 2)
imagesc(P_e_cond_s)
set(gca,'YDir','normal')
colormap jet
xline(find(belt_idx), 'r')
yline(find(belt_idx), 'r')
subplot(2, 2, 1)
plot(mean(P_e_cond_s, 2), 1:length(P_e_cond_s))
subplot(2, 2, 4)
plot(mean(P_e_cond_s, 1))

P_e_cond_s(P_e_cond_s < prctile(P_e_cond_s(:), 95)) = 0;
g = digraph(P_e_cond_s);
% g = digraph(P_se);
LWidths = (g.Edges.Weight).^2;
LWidths = 10 .* LWidths ./ range(LWidths);
LWidths(isnan(LWidths)) = 0;

figure
cm = rbmap('caxis',[0 1]);
% cm = cm(knnsearch(linspace(0, range(g.Edges.Weight), size(cm,1))', g.Edges.Weight), :);
plot(g, 'layout', 'circle', 'LineWidth', 4, 'edgecolor', cm(end, :), 'arrowsize', 0)
axis square


