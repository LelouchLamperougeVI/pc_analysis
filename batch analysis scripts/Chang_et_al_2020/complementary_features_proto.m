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


%% Fig5 a - Reactivated cue/traj features
pc_list = cat(2, rsc.pc_list, ee.pc_list);
num_pc = cellfun(@(x) size(x, 2), clust_stacks{2});

hiepi_psth = cat(1, rsc.hiepi_psth{2}, ee.hiepi_psth{2});
hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = cellfun(@(x) mean(x, 2), hiepi_psth, 'UniformOutput', false);
hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = fast_smooth(hiepi_psth, 5 / 150 * 50);
hiepi_psth = (hiepi_psth - min(hiepi_psth, [], 1)) ./ range(hiepi_psth, 1);

cue_ens = hiepi_psth(:, iscue2 == 1);
traj_ens = hiepi_psth(:, iscue2 == 2);
none_ens = hiepi_psth(:, iscue2 == 0 & num_pc >= 3);

r = corr(cue_ens, traj_ens);
r_null = corr(cue_ens, none_ens);

figure
subplot(1, 3, 1);
[~, idx] = max(cue_ens);
[~, idx] = sort(idx);
imagesc(-cue_ens(:, idx)');
subplot(1, 3, 2);
[~, idx] = min(traj_ens);
[~, idx] = sort(idx);
imagesc(-traj_ens(:, idx)');
subplot(1, 3, 3);
[~, idx] = max(traj_ens);
[~, idx] = sort(idx);
imagesc(-traj_ens(:, idx)');
colormap bone


%% Fig5 b - cue vs traj mean activation strength
y = cat(2, cue_ens, traj_ens);
x = repmat((1:size(cue_ens, 1))', 1, size(y, 2));
g = repelem(1:2, size(x, 1), [size(cue_ens, 2), size(traj_ens, 2)]);

g = gramm('x', x(:), 'y', y(:), 'color', g(:));
g.stat_summary('type', 'sem');
g.axe_property('ylim', [0 .6])
figure
g.draw;


%% Fig5 c - Pearson rho histogram
% g = gramm('x', cat(1, r(:), r_null(:)), 'color', repelem({'traj', 'none'}, [numel(r), numel(r_null)]));
% g.stat_bin('nbins', 50, 'geom', 'stacked_bar', 'normalization', 'probability')
% % g.axe_property('xlim', [-1 1], 'ylim', [0 .05])
% figure
% g.draw;

figure
histogram(r_null(:), linspace(-1, 1, 51), 'Normalization', 'probability', 'DisplayStyle', 'stairs')
hold on
histogram(r(:), linspace(-1, 1, 51), 'Normalization', 'probability', 'DisplayStyle', 'stairs')

signrank(r(:), 0, 'tail', 'left')
skewness(r(:), 0)
[~, p] = kstest(r(:))

signrank(r_null(:), 0, 'tail', 'left')
skewness(r_null(:), 0)
[~, p] = kstest(r_null(:))


%% Reactivated pairs

hiepi_psth = cat(1, rsc.hiepi_psth{2}, ee.hiepi_psth{2});
session = cellfun(@length, hiepi_psth);
session = repelem((1:length(session))', session);

hiepi_psth = cat(2, hiepi_psth{:});
hiepi_psth = cellfun(@(x) mean(x, 2), hiepi_psth, 'UniformOutput', false);
hiepi_psth = cat(2, hiepi_psth{:});
% hiepi_psth = fast_smooth(hiepi_psth, 5 / 150 * 50);
hiepi_psth = (hiepi_psth - min(hiepi_psth, [], 1)) ./ range(hiepi_psth, 1);

hiepi_z = cat(1, rsc.hiepi_z{2}, ee.hiepi_z{2});
temp = {};
for ii = 1:length(hiepi_z)
    temp = cat(1, temp, mat2cell(hiepi_z{ii}', ones(1, size(hiepi_z{ii}, 2))));
end
hiepi_z = temp;

cue_sess = session(iscue2 == 1);
traj_sess = session(iscue2 == 2);
cue_ens = hiepi_psth(:, iscue2 == 1);
traj_ens = hiepi_psth(:, iscue2 == 2);
cue_z = hiepi_z(iscue2 == 1);
traj_z = hiepi_z(iscue2 == 2);

idx = intersect(cue_sess, traj_sess);
% r = nan(length(idx), 1);
ccr = [];
pcc = [];
for ii = 1:length(idx)
%     temp = corr(cue_ens(:, cue_sess == idx(ii)), traj_ens(:, traj_sess == idx(ii)));
%     r(ii) = mean(temp(:));
    for cc = find(cue_sess == idx(ii))'
        for tt = find(traj_sess == idx(ii))'
            t = linspace(0, length(cue_z{cc}) / 19, length(cue_z{cc}));
            a = cue_z{cc};
            b = traj_z{tt};
%             a = fast_smooth(a(:), 5 / 150 * 50);
%             b = fast_smooth(b(:), 5 / 150 * 50);
            temp = isnan(a) | isnan(b);
            a(temp) = [];
            b(temp) = [];
            [r, lags, p] = bxcorr(a, b, t);
            r = r(lags < 1 & lags > -1);
%             [~, pk] = max(r);
%             if r(pk) > p && pk < floor(length(r) /2)
                ccr = cat(2, ccr, r);
                pcc = cat(1, pcc, corr(cue_ens(:, cc), traj_ens(:, tt)));
%             end
        end
    end
end

figure
idx = max(ccr);
[~, idx] = sort(idx);
% imagesc(zscore(ccr(:, idx))');
imagesc(ccr(:, idx)');
xline(19)

figure
boxplot(pcc)
signrank(pcc, 0, 'tail', 'left')


%% try again with onset
hiepi_psth = cat(1, rsc.hiepi_psth{2}, ee.hiepi_psth{2});
session = cellfun(@length, hiepi_psth);
session = repelem((1:length(session))', session);

hiepi_ens = cat(2, hiepi_psth{:});
hiepi_ens = cellfun(@(x) mean(x, 2), hiepi_ens, 'UniformOutput', false);
hiepi_ens = cat(2, hiepi_ens{:});
hiepi_ens = fast_smooth(hiepi_ens, 5 / 150 * 50);
hiepi_ens = (hiepi_ens - min(hiepi_ens, [], 1)) ./ range(hiepi_ens, 1);

hiepi_psth = hiepi_psth(~cellfun(@isempty, hiepi_psth));
hiepi_psth = cat(2, hiepi_psth{:});
cue_psth = hiepi_psth(iscue2 == 1);
traj_psth = hiepi_psth(iscue2 == 2);

cue_stacks = clust_stacks{2}(iscue2 == 1);
traj_stacks = clust_stacks{2}(iscue2 == 2);

hiepi_on = cat(1, rsc.hiepi_lfp_pw{2}, ee.hiepi_lfp_pw{2});
temp = {};
for ii = 1:length(hiepi_on)
    for jj = 1:length(hiepi_on{ii})
        temp = cat(1, temp, hiepi_on{ii}{jj}.on);
    end
end
hiepi_on = temp;
cue_on = hiepi_on(iscue2 == 1);
traj_on = hiepi_on(iscue2 == 2);

cue_sess = session(iscue2 == 1);
traj_sess = session(iscue2 == 2);
cue_ens = hiepi_ens(:, iscue2 == 1);
traj_ens = hiepi_ens(:, iscue2 == 2);
cue_z = hiepi_z(iscue2 == 1);
traj_z = hiepi_z(iscue2 == 2);

idx = intersect(cue_sess, traj_sess);

delta = [];
pcc = [];

count = 1;
master_idx = [];
for ii = 1:length(idx)
    for cc = find(cue_sess == idx(ii))'
        for tt = find(traj_sess == idx(ii))'
            d = cue_on{cc}(:) - traj_on{tt}(:)';
            if size(d, 1) < size(d, 2)
                d = d';
            end
            [~, temp] = min(abs(d));
            temp = d(sub2ind(size(d), temp, 1:size(d, 2)));
            temp = sum(abs(temp) < 0.5) / length(temp);
%             temp = median(temp);
%             temp(abs(temp) > 1) = [];
%             [~, temp2] = sort(abs(temp));
%             temp = mean(temp(temp2(1:round(length(temp2)*.2))));
            delta = cat(1, delta, temp);
            
            pcc = cat(1, pcc, corr(cue_ens(:, cc), traj_ens(:, tt)));
            
            master_idx = cat(1, master_idx, [cc, tt]);
            
%             if delta(end) < .5
%                 if ~mod(count - 1, 25)
%                     figure
%                 end
%                 subplot(5, 5, mod(count - 1, 25) + 1);
%                 plot(cue_ens(:, cc))
%                 hold on
%                 plot(traj_ens(:, tt))
%                 title(num2str(pcc(end)))
%                 count = count + 1;
%             end
        end
    end
end

pcc(isnan(delta)) = [];
delta(isnan(delta)) = [];

temp = cellfun(@(x) size(x, 2), traj_stacks(master_idx(:, 2)));
% idx = find(pcc < -.4 & delta > .1);
idx = find(pcc < -.4 & temp > 6);
count = 1;
for ii = 1:length(idx)
    if ~mod(count - 1, 5)
        figure
    end
    subplot(3, 5, mod(count - 1, 5) + 1);
    stack = cue_stacks{master_idx(idx(ii), 1)};
    [~, order] = max(stack);
    [~, order] = sort(order);
    imagesc(stack(:, order)');
    
    subplot(3, 5, mod(count - 1, 5) + 1 + 5);
    stack = traj_stacks{master_idx(idx(ii), 2)};
    [~, order] = max(stack);
    [~, order] = sort(order);
    imagesc(stack(:, order)');
    
    subplot(3, 5, mod(count - 1, 5) + 1 + 10);
    plot(cue_ens(:, master_idx(idx(ii), 1)))
    hold on
    plot(traj_ens(:, master_idx(idx(ii), 2)))
    count = count + 1;
    xlim([1 50])
end



