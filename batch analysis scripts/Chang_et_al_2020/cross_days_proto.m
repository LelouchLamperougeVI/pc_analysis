%% classify them ensembles... again...

clear all
load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clust_stacks', 'whole_stack', 'clusts');

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

iscue1 = mat2cell(iscue1, cellfun(@length, clusts{1}));
iscue2 = mat2cell(iscue2, cellfun(@length, clusts{2}));

clearvars -except iscue1 iscue2 whole_stack clusts

%%

% clear all

% iscue2 = load('/mnt/storage/HaoRan/RRR_motor/M2/manual_class.mat');
% iscue2 = iscue2.temp;

% load('/mnt/storage/HaoRan/RRR_motor/M2/auto_class.mat');

% load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clusts', 'whole_stack')

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
% animals = dir(fullfile(root, 'RSC*'));

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals([animals.isdir]).name};

anal_f = cell(length(animals), 1); % grosse cochonne, tu pensais a quoi? ;)
daydate = cell(length(animals), 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    anal_f{a} = cellfun(@(x) fullfile(root, animals{a}, x, 'analysis.mat'), sessions, 'uniformoutput', false);
    
    daydate{a} = cellfun(@(x) datetime(x, 'InputFormat', 'yyyy_MM_dd'), sessions);
end

idx = cellfun(@length, daydate) > 2;
for ii = find(idx(:)')
    daydate{ii} = [0 cumsum(split((caldiff(daydate{ii})), {'days'}))];
end
for ii = find(~idx(:)')
    if ~isempty(daydate{ii})
        daydate{ii} = 0;
    else
        daydate{ii} = [];
    end
    anal_f{ii} = '';
end

temp = daydate(cellfun(@length, daydate) > 1);
x = cellfun(@(x) 1:length(x), temp, 'UniformOutput', false);
y = temp;
g = 1:length(temp);
cmap = distinguishable_colors(length(temp), 'w');
g = gramm('x', x, 'y', y, 'color', g);
g.set_color_options('map', cmap);
g.geom_line('alpha', .25, 'dodge', .5);
g.geom_point('alpha', .4, 'dodge', .5);
g.set_names('x', 'consecutive recording number', 'y', 'days passed since first recording', 'color', 'animal');
g.axe_property('XLim', [0 11], 'ylim', [-1 41]);
figure
g.draw();

% figure
% [c, e] = histcounts(cell2mat(daydate'));
% bar(e(1:end-1)+.5, c)


%% Consecutive days mode (this block is mutually exclusive with next block)
target_days = [0 1];

idx = cell(length(daydate), 1);
for ii = 1:length(daydate)
    idx{ii} = [(1:(length(daydate{ii})-1))', (2:length(daydate{ii}))'];
    anal_f{ii} = anal_f{ii}(idx{ii});
    temp = iscue2(sum(cellfun(@length, daydate(1:(ii-1))))+1 : sum(cellfun(@length, daydate(1:ii)))-1)';
    anal_f{ii}(~cellfun(@any, temp), :) = [];
    idx{ii}(~cellfun(@any, temp), :) = [];
    idx{ii} = idx{ii} + sum(cellfun(@length, daydate(1:(ii-1))));
end

idx = cell2mat(idx);
anal_f = cat(1, anal_f{:});


%% Choose targets for cross days (fixed day delta mode)
target_days = [0 2];

idx = cellfun(@(x) all(ismember(target_days, x)), daydate);
temp = [];
for ii = 1:length(daydate)
    if idx(ii)
        temp = cat(1, temp, sum((daydate{ii}(:) == target_days) .* (1:length(target_days)), 2));
    else
        temp = cat(1, temp, zeros(length(daydate{ii}), 1));
    end
end

daydate = temp;

idx = daydate(:) == 1:length(target_days);
idx = arrayfun(@(x) find(idx(:, x)), 1:length(target_days), 'UniformOutput', false);
idx = cell2mat(idx);


%% Register them masks and find overlapping ROIs
% Quand la realite frappe a la porte... t'as rien a faire que de te
% preparer a etre insere :)

olap_thres = .5;

load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'masks');
masks = masks{2};

ROIs = cell(length(target_days), length(target_days), size(idx, 1)); % rows x cols x session

% this code was stolen from multiplane/register.m
for s = 1:size(idx, 1)
    for r = 1:length(target_days)
        for c = 1:length(target_days)
            maskReg = masks{idx(s, r)};
            if r == c
                ROIs{r, c, s} = 1:max(maskReg(:));
                continue
            end
            
            maskA = regmasks(maskReg, masks{idx(s, c)});
            
            stack = maskReg == permute(1:max(maskReg(:)), [1 3 2]); % tried to be lazy, ended up overcomplicating things...
            stack = maskA .* stack;
            stack = reshape(stack, [numel(maskA) size(stack,3)]); % le gros machin complique a transformations des matrices irregulieres
            stack = arrayfun(@(x) accumarray(stack(:, x)+1, 1, [max(maskA(:))+1 1]), 1:size(stack,2), 'uniformoutput',false);
            stack = cell2mat(stack);
            stack = stack(2:end,:);
            
            movPix = accumarray(maskReg(:)+1, 1, [max(maskReg(:))+1 1]);
            movPix = movPix(2:end);
            
            fixPix = accumarray(maskA(:)+1, 1, [max(maskA(:))+1 1]); % number of pixels for fixed ROIs
            fixPix = fixPix(2:end);
            
            stack = stack ./ ( fixPix + movPix' - stack );
            roi = arrayfun(@(x) find(stack(:,x) > olap_thres), 1:size(stack,2), 'uniformoutput',false);
            
            roi(cellfun(@isempty, roi)) = {nan};
            roi = cell2mat(roi);
            
            ROIs{r, c, s} = roi;
        end
    end
end
sess_idx = idx;


%% make cross days stacks

% putain de merde... ce n'etait pas si difficile la premiere fois que je
% l'aie realise pourtant...

cross_stack = cell(length(target_days), length(target_days));
for r = 1:length(target_days)
    for c = 1:length(target_days)
        temp = [];
        for s = 1:size(sess_idx, 1)
            idx = ROIs{r, c, s}(cell2mat(clusts{2}{sess_idx(s, r)}(iscue2{sess_idx(s, r)}==1)));
            temptemp = nan(size(whole_stack{sess_idx(s, c)}, 1), length(idx));
            temptemp(:, ~isnan(idx)) = whole_stack{sess_idx(s, c)}(:, idx(~isnan(idx)));
            temp = cat(2, temp, temptemp);
        end
        cross_stack{r, c} = temp;
    end
end

figure
for ii = 1:length(target_days)
    [~, order] = max(cross_stack{ii, ii});
    [~, order] = sort(order);
    idx = cellfun(@(x) any(isnan(x)), cross_stack(ii, :), 'UniformOutput', false);
    idx = find(any(cell2mat(idx')));
    order(ismember(order, idx)) = [];
    disp(numel(order))
    for jj = 1:length(target_days)
        subplot(length(target_days), length(target_days), (ii-1)*length(target_days) + jj)
        imagesc(cross_stack{ii, jj}(:, order)')
        colormap viridis
    end
end

figure
r = corr(cross_stack{1, 1}(:, ~any(isnan(cross_stack{1, 2})))', cross_stack{1, 2}(:, ~any(isnan(cross_stack{1, 2})))');
imagesc(r)
axis square
caxis([-1 1])
colorbar
rbmap('interp', 257);

r_cue = diag(corr(cross_stack{1, 1}, cross_stack{1, 2}));


cross_stack = cell(length(target_days), length(target_days));
for r = 1:length(target_days)
    for c = 1:length(target_days)
        temp = [];
        for s = 1:size(sess_idx, 1)
            idx = ROIs{r, c, s}(cell2mat(clusts{2}{sess_idx(s, r)}(iscue2{sess_idx(s, r)}==2)));
            temptemp = nan(size(whole_stack{sess_idx(s, c)}, 1), length(idx));
            temptemp(:, ~isnan(idx)) = whole_stack{sess_idx(s, c)}(:, idx(~isnan(idx)));
            temp = cat(2, temp, temptemp);
        end
        cross_stack{r, c} = temp;
    end
end

figure
for ii = 1:length(target_days)
    [~, order] = max(cross_stack{ii, ii});
    [~, order] = sort(order);
    idx = cellfun(@(x) any(isnan(x)), cross_stack(ii, :), 'UniformOutput', false);
    idx = find(any(cell2mat(idx')));
    order(ismember(order, idx)) = [];
    disp(numel(order))
    for jj = 1:length(target_days)
        subplot(length(target_days), length(target_days), (ii-1)*length(target_days) + jj)
        imagesc(cross_stack{ii, jj}(:, order)')
        colormap viridis
    end
end

figure
r = corr(cross_stack{1, 1}(:, ~any(isnan(cross_stack{1, 2})))', cross_stack{1, 2}(:, ~any(isnan(cross_stack{1, 2})))');
imagesc(r)
axis square
caxis([-1 1])
colorbar
rbmap('interp', 257);

r_traj = diag(corr(cross_stack{1, 1}, cross_stack{1, 2}));

% sem = @(x) std(x, 'omitnan') / sum(~isnan(x));
% sem = @(x) std(x, 'omitnan');
% 
% figure
% bar([mean(r_cue, 'omitnan'), mean(r_traj, 'omitnan')], 'k');
% hold on
% errorbar([mean(r_cue, 'omitnan'), mean(r_traj, 'omitnan')], [sem(r_cue), sem(r_traj)], 'k', 'LineStyle', 'none');
% [~, pval] = ttest2(r_cue, r_traj);
% title(['two-sample two-tailed t-test; p = ' num2str(pval)])

figure
boxplot([r_cue; r_traj], repelem((1:2)', [numel(r_cue); numel(r_traj)], 1));
p = ranksum(r_cue, r_traj, 'tail', 'right');
title(['two-sample one-tailed ranksum; p = ' num2str(p)])




%% Quantify silent cells
% L'intention ici est de voir si les cellules qui représentent des
% trajectoires deviennent silencieuses dans le prochain jour d'acquisition.
% Cette hypothèse est basée sur la notion dont les neurones pyramidaux de
% l'hippocampe réorganisent plus couramment et se privent de l'opportunité
% de s'attacher aux ressources externes, sous peine de perdre la capacité
% à encoder de nouvelles informations qu'ils reçoivent en ligne.

silent_cues = cell(1, size(sess_idx, 1));
silent_traj = cell(1, size(sess_idx, 1));
silent_all = cell(1, size(sess_idx, 1));

for s = 1:size(sess_idx, 1)
    idx = clusts{2}{sess_idx(s, 1)}(iscue2{sess_idx(s, 1)}==1); % la meme chose que plutot....
    silent_cues{s} = cellfun(@(x) sum(isnan(ROIs{1, 2, s}(x))) / numel(x), idx);
    idx = clusts{2}{sess_idx(s, 1)}(iscue2{sess_idx(s, 1)}==2);
    silent_traj{s} = cellfun(@(x) sum(isnan(ROIs{1, 2, s}(x))) / numel(x), idx);
    silent_all{s} = sum(isnan(ROIs{1, 2, s})) / numel(ROIs{1, 2, s});
end

silent_cues = arrayfun(@(x) silent_cues{x} ./ silent_all{x}, 1:length(silent_all), 'UniformOutput', false);
silent_traj = arrayfun(@(x) silent_traj{x} ./ silent_all{x}, 1:length(silent_all), 'UniformOutput', false);

silent_cues = cell2mat(silent_cues);
silent_traj = cell2mat(silent_traj);

figure
cdfplot(silent_cues);
hold on
cdfplot(silent_traj);
% Il semblerait que la fraction de superposition varie en fonction du jour,
% tandis qu'elle est nullement influencée par l'identité de l'ensemble.


%% Population vector stability
shuffles = 1e3;

cross_stack = cell(length(target_days), length(target_days));
for r = 1:length(target_days)
    for c = 1:length(target_days)
        temp = [];
        for s = 1:size(sess_idx, 1)
            idx = ROIs{r, c, s}(cell2mat(clusts{2}{sess_idx(s, r)}(iscue2{sess_idx(s, r)}==1)));
            temptemp = nan(size(whole_stack{sess_idx(s, c)}, 1), length(idx));
            temptemp(:, ~isnan(idx)) = whole_stack{sess_idx(s, c)}(:, idx(~isnan(idx)));
            temp = cat(2, temp, temptemp);
        end
        cross_stack{r, c} = temp;
    end
end

cross_stack = cross_stack(1, :);
idx = find(~any(isnan(cross_stack{2})));
r = corr(cross_stack{1}(:, idx)', cross_stack{2}(:, idx)');
r_cue = diag(r);

r_cue_boot = zeros(shuffles, length(r_cue));
for ii = 1:shuffles
    boot = randsample(idx, length(idx), true);
    r = corr(cross_stack{1}(:, boot)', cross_stack{2}(:, boot)');
    r_cue_boot(ii, :) = diag(r);
end

cross_stack = cell(length(target_days), length(target_days));
for r = 1:length(target_days)
    for c = 1:length(target_days)
        temp = [];
        for s = 1:size(sess_idx, 1)
            idx = ROIs{r, c, s}(cell2mat(clusts{2}{sess_idx(s, r)}(iscue2{sess_idx(s, r)}==2)));
            temptemp = nan(size(whole_stack{sess_idx(s, c)}, 1), length(idx));
            temptemp(:, ~isnan(idx)) = whole_stack{sess_idx(s, c)}(:, idx(~isnan(idx)));
            temp = cat(2, temp, temptemp);
        end
        cross_stack{r, c} = temp;
    end
end

cross_stack = cross_stack(1, :);
idx = find(~any(isnan(cross_stack{2})));
r = corr(cross_stack{1}(:, idx)', cross_stack{2}(:, idx)');
r_traj = diag(r);

r_traj_boot = zeros(shuffles, length(r_cue));
for ii = 1:shuffles
    boot = randsample(idx, length(idx), true);
    r = corr(cross_stack{1}(:, boot)', cross_stack{2}(:, boot)');
    r_traj_boot(ii, :) = diag(r);
end


ci = cat(2, r_cue_boot, r_traj_boot);
ci = prctile(ci, [2.5 97.5]);
r = cat(2, r_cue(:)', r_traj(:)');
x = linspace(0, 150, length(r_cue));
x = [x, x];
g = repelem({'cue', 'traj'}, [length(r_cue), length(r_traj)]);
g = gramm('x', x, 'y', r, 'color', g, 'ymin', ci(1, :), 'ymax', ci(2, :));
g.geom_line;
g.geom_interval;
g.axe_property('YLim', [0 1]);
figure
g.draw;



%% Cross days bayesian decoding fun :)

e = [];
n = [];
for s = 1:size(sess_idx, 1)
    a1 = load(anal_f{s, 1}); a1 = a1.analysis;
    a2 = load(anal_f{s, 2}); a2 = a2.analysis;
    
    idx = ROIs{1, 2, s};
%     idx( setxor(1:length(idx), cell2mat(clusts{2}{sess_idx(s, 1)}(iscue2{sess_idx(s, 1)} == 2)) )) = nan;
    
    if all(isnan(idx))
        continue
    end

    [decoded, P, pos, err] = bayes_paired(a1, a2, idx);
    e = cat(1, e, err(:, 1)');
    n = cat(1, n, sum(~isnan(idx)));
end


[decoded, P, pos, err] = bayes_paired(a1, a1, 1:size(a1.deconv, 2)); % troubleshooting...

    idx = clusts{2}{sess_idx(s, 1)}(iscue2{sess_idx(s, 1)}==1); % la meme chose que plutot....
    silent_cues{s} = cellfun(@(x) sum(isnan(ROIs{1, 2, s}(x))) / numel(x), idx);
    idx = clusts{2}{sess_idx(s, 1)}(iscue2{sess_idx(s, 1)}==2);
    silent_traj{s} = cellfun(@(x) sum(isnan(ROIs{1, 2, s}(x))) / numel(x), idx);
    silent_all{s} = sum(isnan(ROIs{1, 2, s})) / numel(ROIs{1, 2, s});
% end




