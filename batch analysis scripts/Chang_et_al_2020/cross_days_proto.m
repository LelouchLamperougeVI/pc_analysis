clear all

iscue2 = load('/mnt/storage/HaoRan/RRR_motor/M2/manual_class.mat');
iscue2 = iscue2.temp;

% load('/mnt/storage/HaoRan/RRR_motor/M2/auto_class.mat');

load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'clusts', 'whole_stack')

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
% animals = dir(fullfile(root, 'RSC*'));

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

daydate = cell(length(animals), 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
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
end

figure
hold on
for ii = 1:length(daydate); plot(daydate{ii}); end
xlabel('next recording day');
ylabel('days passed since first recording')

figure
[c, e] = histcounts(cell2mat(daydate'));
bar(e(1:end-1)+.5, c)


%% Consecutive days mode (this block is mutually exclusive with next block)
target_days = [0 1];

idx = cell(length(daydate), 1);
for ii = 1:length(daydate)
    idx{ii} = [(1:(length(daydate{ii})-1))', (2:length(daydate{ii}))'];
    temp = iscue2(sum(cellfun(@length, daydate(1:(ii-1))))+1 : sum(cellfun(@length, daydate(1:ii)))-1)';
    idx{ii}(~cellfun(@any, temp), :) = [];
    idx{ii} = idx{ii} + sum(cellfun(@length, daydate(1:(ii-1))));
end

idx = cell2mat(idx);


%% Choose targets for cross days
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
olap_thres = .5;

load('/mnt/storage/rrr_magnum/M2/hiepi4.mat', 'masks');
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
            stack = reshape(stack, [numel(maskA) size(stack,3)]);
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
    for jj = 1:length(target_days)
        subplot(length(target_days), length(target_days), (ii-1)*length(target_days) + jj)
        imagesc(-cross_stack{ii, jj}(:, order)')
        colormap gray
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
    for jj = 1:length(target_days)
        subplot(length(target_days), length(target_days), (ii-1)*length(target_days) + jj)
        imagesc(-cross_stack{ii, jj}(:, order)')
        colormap gray
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

sem = @(x) std(x, 'omitnan') / sum(~isnan(x));
sem = @(x) std(x, 'omitnan');

figure
bar([mean(r_cue, 'omitnan'), mean(r_traj, 'omitnan')], 'k');
hold on
errorbar([mean(r_cue, 'omitnan'), mean(r_traj, 'omitnan')], [sem(r_cue), sem(r_traj)], 'k', 'LineStyle', 'none');
[~, pval] = ttest2(r_cue, r_traj);
title(['two-sample two-tailed t-test; p = ' num2str(pval)])


