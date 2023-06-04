clear all

mf = 60;
tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro');

root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));

% root = '/mnt/storage/rrr_magnum/M2';
% animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

count = 1;

results = []; %EV, REV, e1, e2
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
        rest1.hPICA;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
        rest2.hPICA;
        
        du = rest1.analysis.original_deconv; % run
        d1 = rest1.twop.deconv;
        d1(any(isnan(d1), 2), :) = [];
        d2 = rest2.twop.deconv;
        d2(any(isnan(d2), 2), :) = [];

        e1 = nan(mf, length(rest1.ensembles.clust));
        e2 = nan(mf, length(rest2.ensembles.clust));
        EV = nan(mf, 1);
        REV = nan(mf, 1);
        for cf = 1:mf
            dt = downsample(movmean(du, cf), cf);
%             dc1 = downsample(movmean(d1, cf), cf);
%             dc2 = downsample(movmean(d2, cf), cf);
            dc1 = d1;
            dc2 = d2;

%             [EV(cf), REV(cf), e1(:, :, cf), e2(:, :, cf)] = ev(rest1, d, rest2);

            dt = (dt - mean(dt)) ./ std(dt); % z-score X to null mean and unit variance
            dc1 = (dc1 - mean(dc1)) ./ std(dc1); % z-score X to null mean and unit variance
            dc2 = (dc2 - mean(dc2)) ./ std(dc2); % z-score X to null mean and unit variance
            
            try
                EV(cf) = tev(dt, rest2.hiepi.pc) / tev(dc2, rest2.hiepi.pc);
%                 EV(cf) = tev(dt, rest2.hiepi.pc);
            catch
            end
            try
                REV(cf) = tev(dt, rest1.hiepi.pc) / tev(dc1, rest2.hiepi.pc);
%                 REV(cf) = tev(dt, rest1.hiepi.pc);
            catch
            end
            try
                e1(cf, :) = arrayfun(@(ii) tev(dt, rest1.hiepi.pc(:, ii)), 1:size(rest1.hiepi.pc, 2)) ./ arrayfun(@(ii) tev(dc1, rest1.hiepi.pc(:, ii)), 1:size(rest1.hiepi.pc, 2));
%                 e1(cf, :) = arrayfun(@(ii) tev(dt, rest1.hiepi.pc(:, ii)), 1:size(rest1.hiepi.pc, 2));
            catch
            end
            try
                e2(cf, :) = arrayfun(@(ii) tev(dt, rest2.hiepi.pc(:, ii)), 1:size(rest2.hiepi.pc, 2)) ./ arrayfun(@(ii) tev(dc2, rest2.hiepi.pc(:, ii)), 1:size(rest2.hiepi.pc, 2));
%                 e2(cf, :) = arrayfun(@(ii) tev(dt, rest2.hiepi.pc(:, ii)), 1:size(rest2.hiepi.pc, 2));
            catch
            end
        end
%         
%         e1 = permute(e1(:, 2, :), [3, 1, 2]);
%         e2 = permute(e2(:, 1, :), [3, 1, 2]);
%         
        results = cat(1, results, {EV, REV, e1, e2});
        
        count = count + 1;
    end
end

%%
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'clust_stacks', 'clusts', 'whole_stack');

% classify ensembles
clust_stacks{1} = cat(1, rsc.clust_stacks{1}, ee.clust_stacks{1});
clust_stacks{2} = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});
% clust_stacks{1} = rsc.clust_stacks{1};
% clust_stacks{2} = rsc.clust_stacks{2};

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

clearvars -except iscue1 iscue2


%%
results1 = load('/mnt/storage/HaoRan/RRR_motor/M2/review_time_compression.mat');
results2 = load('/mnt/storage/rrr_magnum/M2/review_time_compression.mat');
EV = cat(2, results1.results{:, 1}, results2.results{:, 1});
e2 = cat(2, results1.results{:, 3}, results2.results{:, 3}, results1.results{:, 4}, results2.results{:, 4});
% e2 = cat(2, results1.results{:, 4}, results2.results{:, 4});
e2 = -sqrt((e2 - 1) .^ 2) + 1;

% iscue = iscue2;
iscue = cat(1, iscue1, iscue2);

figure
% subplot(1, 3, 1)
% errorbar(mean(zscore(EV), 2), sem(zscore(EV), 2));
% xlabel('compression factor')
% ylabel('explained variance');
% title('session total');

subplot(2, 3, 1:3)
% errorbar(mean(zscore(e2(:, iscue == 1)), 2), sem(zscore(e2(:, iscue == 1)), 2));
errorbar(mean(e2(:, iscue == 1), 2), sem(e2(:, iscue == 1), 2), 'Color', hex2dec(["00" "8b" "52"]) ./ 255);
hold on
% yyaxis right
% errorbar(mean(zscore(e2(:, iscue == 2)), 2), sem(zscore(e2(:, iscue == 2)), 2));
errorbar(mean(e2(:, iscue == 2), 2), sem(e2(:, iscue == 2), 2), 'Color', hex2dec(["41" "25" "54"]) ./ 255);
xlabel('compression factor')
ylabel('divergence');
legend('cue', 'traj');

subplot(2, 3, 4:5);
[~, pks] = max(e2, [], 1);
h = cdfplot(pks(iscue == 1));
h.Color = hex2dec(["00" "8b" "52"]) ./ 255;
hold on
h = cdfplot(pks(iscue == 2));
h.Color = hex2dec(["41" "25" "54"]) ./ 255;
legend('cue', 'traj');
xlabel('optimal cf');
ylabel('cum. freq.');

subplot(2, 3, 6);
boxplot(pks(~(iscue == 0)), iscue(~(iscue == 0)), 'PlotStyle', 'compact',...
    'Colors', cat(1, hex2dec(["00" "8b" "52"]) ./ 255, hex2dec(["41" "25" "54"]) ./ 255), ...
    'Labels', {'cue', 'traj'});
ylabel('optimal cf')

[~, p] = kstest2(pks(iscue == 1), pks(iscue == 2))
p = ranksum(pks(iscue == 1), pks(iscue == 2))



%%
[W, H] = nnmf(e2, 3);

H = H(:, ~(iscue == 0))';
iscue_ = iscue(~(iscue == 0));
g = kmeans(H, 2);

if sum(g == iscue_) < (length(g) / 2)
    g = -(g - 1) + 2;
end

h = figure;
ax(1) = subplot(1, 3, 1);
plot(W);
title('NNMF features');
xlabel('compression factor');
ylabel('factorised divergence');
legend('feat 1', 'feat 2', 'feat 3');

ax(1) = subplot(1, 3, 2);
scatter3(H(iscue_ == 1, 1), H(iscue_ == 1, 2), H(iscue_ == 1, 3), 'MarkerEdgeColor', hex2dec(["00" "8b" "52"]) ./ 255, 'Marker', 'o')
hold on
scatter3(H(iscue_ == 2, 1), H(iscue_ == 2, 2), H(iscue_ == 2, 3), 'MarkerEdgeColor', hex2dec(["41" "25" "54"]) ./ 255, 'Marker', 'o')
title('criteria labeled');
xlabel('feat 1');
ylabel('feat 2');
zlabel('feat 3');
legend('cue', 'traj');

ax(2) = subplot(1, 3, 3);
scatter3(H(g == 1, 1), H(g == 1, 2), H(g == 1, 3), 'MarkerEdgeColor', hex2dec(["00" "8b" "52"]) ./ 255, 'Marker', 'o')
hold on
scatter3(H(g == 2, 1), H(g == 2, 2), H(g == 2, 3), 'MarkerEdgeColor', hex2dec(["41" "25" "54"]) ./ 255, 'Marker', 'o')
title('k-means clustered');
xlabel('feat 1');
ylabel('feat 2');
zlabel('feat 3');
legend('cue', 'traj');

Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(h, 'StoreTheLink', Link);

figure
confusionchart(accumarray([iscue_, g], 1), {'cue', 'traj'}, 'Title', 'k-means clustering', ...
    'ColumnSummary', 'column-normalized', 'RowSummary', 'row-normalized', 'Normalization', 'total-normalized')

chisq2(g, iscue_)
