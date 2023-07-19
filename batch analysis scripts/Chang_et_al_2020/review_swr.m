clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'hiepi_struct', 'pc_list', 'file', 'swr_struct');

% classify ensembles
clust_stacks = ee.clust_stacks;
% clust_stacks{2}([13:49, 223:256]) = []; %%%%%%%%%%%%%%%%

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

clearvars -except ee iscue1 iscue2

% params
wdw = 1;
bandwidth = .1;

% idx = double(cellfun(@isempty, ee.clusts{2})); %%%%%%%%%%%%%%%%
% idx([5:7, 41:44]) = 2;
% idx(idx == 1) = [];
% ee.hiepi_struct{2}(idx == 2) = [];
% 
% ee.swr_struct{2}([5:7, 41:44]) = []; %%%%%%%%%%%%%%%%
% ee.clusts{2}([5:7, 41:44]) = [];
% ee.pc_list([5:7, 41:44]) = [];

rest = 2;
swr_on = {ee.swr_struct{rest}(:).swr_on};
swr_on = swr_on(:);

idx = cellfun(@isempty, ee.clusts{rest});
hiepi_on = cell(length(idx), 1);
hiepi_on(~idx) = ee.hiepi_struct{rest}(:);

xcf = [];
for ii = 1:length(swr_on)
    for jj = 1:length(hiepi_on{ii})
        temp = hiepi_on{ii}(jj).on(:)' - swr_on{ii}(:);
        xcf = cat(1, xcf, {temp(:)});
    end
end

xi = linspace(-wdw, wdw, 1e3);
f = cellfun(@(x) ksdensity(x, xi, 'bandwidth', bandwidth), xcf, 'UniformOutput', false);
f = cell2mat(f);
[~, delay] = max(f, [], 2);
delay = xi(delay);


load('/mnt/storage/rrr_magnum/M2/ccpc_ll.mat', 'md2');
% md2([5:7, 41:44]) = []; %%%%%%%%%%%%%%%%

clusts = ee.clusts{2};
count = 1;
for ii = 1:length(clusts)
    for jj = 1:length(clusts{ii})
        temp = md2{ii}.ll(intersect(clusts{ii}{jj}, ee.pc_list{ii}), :);
%         temp = md2{ii}.ll(clusts{ii}{jj}, :);
        ratio{count} = -2 .* (temp(:, 2) - temp(:, 1));
%         ssize{count} = md{ii}.ssize(clusts{ii}{jj});
        count = count + 1;
    end
end

[~, pks] = max(f, [], 2);
% ll = cellfun(@mean, ratio)';
% ll(cellfun(@length, ratio) < 3) = nan;
ll = cat(1, ratio{:});
pks = xi(pks);
pks = repelem(pks, cellfun(@length, ratio));


figure
[~, order] = max(f, [], 2);
[~, order] = sort(order);

subplot(2, 2, [1, 3]);
imagesc(zscore(f(order, :)')')
colormap jet

subplot(2, 2, 2);
histogram(pks, linspace(-1, 1, 26), 'Normalization', 'probability')

x = ll;
y = pks;
xx = linspace(min(x), max(x), 1e4);

[b, stats] = robustfit(x, y);

n = length(x);
fit = b(1) + x .* b(2);
tcrit = tinv([.025, .975], n - 2);
se = sqrt( sum((y(:) - fit(:) ).^2) / (n - 2) ) .* ...
    sqrt(1/n + (xx - mean(x)).^2 ./ sum((x - mean(x)).^2));

rsqr = 1 - sum((y(:) - fit(:)).^2) / sum((y(:) - mean(y)).^2);
rsqr_adj = 1 - (1 - rsqr) * (n - 1) / (n - length(b));

fit = b(1) + xx(:) .* b(2);
ci = fit + se(:) .* tcrit(:)';

subplot(2, 2, 4);
% scatter(x, y, 'MarkerFaceAlpha', .1, 'MarkerFaceColor', 'r', 'MarkerEdgeAlpha', 0);
hold on
plot(xx, ci, 'k--')
plot(xx, fit, 'r');

title({'robust linear regression', ['y ~ ' num2str(b(1)) ' + x * ' num2str(b(2))], ...
    ['b0 t=' num2str(stats.t(1)) '; p=' num2str(stats.p(1))], ...
    ['b1 t=' num2str(stats.t(2)) '; p=' num2str(stats.p(2))], ...
    ['residuals df = ' num2str(stats.dfe)]});
xlabel('likelihood ratio');
ylabel('time from SWR (sec)');



