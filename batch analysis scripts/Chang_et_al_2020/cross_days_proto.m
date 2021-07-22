clear all

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

% figure
% hold on
% for ii = 1:length(daydate); plot(daydate{ii}); end
% xlabel('next recording day');
% ylabel('days passed since first recording')
%
% figure
% [c, e] = histcounts(cell2mat(daydate'));
% bar(e(1:end-1)+.5, c)

%% Choose targets for cross days
target_days = [0 2 8];

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


%% Register them masks and find overlapping ROIs
olap_thres = .5;

load('/mnt/storage/rrr_magnum/M2/hiepi4.mat', 'masks');
masks = masks{2};

idx = daydate(:) == 1:length(target_days);
idx = arrayfun(@(x) find(idx(:, x)), 1:length(target_days), 'UniformOutput', false);
idx = cell2mat(idx);

ROIs = cell(length(target_days), length(target_days), size(idx, 1)); % rows x cols x session

% that code was stolen from multiplane/register.m
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




