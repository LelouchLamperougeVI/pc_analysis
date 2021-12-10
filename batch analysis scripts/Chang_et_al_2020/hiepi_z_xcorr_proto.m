clear all

load('/mnt/storage/rrr_magnum/M2/hiepi4.mat');
% load('/mnt/storage/HaoRan/RRR_motor/M2/hiepi4.mat');

fr_thres = .5;
traj_thres = 3; %min number of pc per ensemble
l_thres = 30; % length to be considered cue ensemble

l1 = cell(length(swr_clust_stack{1}), 1); s1=l1; e1=l1;
for s = 1:length(swr_clust_stack{1})
    l1{s} = cell(length(swr_clust_stack{1}{s}), 1);
    for c = 1:length(swr_clust_stack{1}{s})
        stack = swr_clust_stack{1}{s}{c};
        traj = any(stack > fr_thres, 2);
        [~, starts, ends] = traj_length(traj, 1);

        stack = repmat(stack, 2, 1);
        idx = false(length(starts{1}), 1);
        for t = 1:length(starts{1})
            temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
            temp = any(temp > fr_thres, 1);
            idx(t) = sum(temp) < traj_thres;
        end

        [l1{s}{c}, s1{s}{c}, e1{s}{c}] = traj_length(traj);
        l1{s}{c} = l1{s}{c}{1}(~idx);
        s1{s}{c} = s1{s}{c}{1}(~idx);
        e1{s}{c} = e1{s}{c}{1}(~idx);
    end
end

l2 = cell(length(swr_clust_stack{2}), 1); s2=l2; e2=l2;
for s = 1:length(swr_clust_stack{2})
    l2{s} = cell(length(swr_clust_stack{2}{s}), 1);
    for c = 1:length(swr_clust_stack{2}{s})
        stack = swr_clust_stack{2}{s}{c};
        traj = any(stack > fr_thres, 2);
        [~, starts, ends] = traj_length(traj, 1);

        stack = repmat(stack, 2, 1);
        idx = false(length(starts{1}), 1);
        for t = 1:length(starts{1})
            temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
            temp = any(temp > fr_thres, 1);
            idx(t) = sum(temp) < traj_thres;
        end

        [l2{s}{c}, s2{s}{c}, e2{s}{c}] = traj_length(traj);
        l2{s}{c} = l2{s}{c}{1}(~idx);
        s2{s}{c} = s2{s}{c}{1}(~idx);
        e2{s}{c} = e2{s}{c}{1}(~idx);
    end
end

belt;
iscue1 = cell(length(l1), 1);
for ii = 1:length(l1) %classify cue/traj ensembles
    iscue1{ii} = zeros(length(l1{ii}), 1);
    for jj = 1:length(l1{ii})
        if isempty(l1{ii}{jj}); continue; end
        temp = any(s1{ii}{jj} < cue_centres & e1{ii}{jj} > cue_centres & l1{ii}{jj} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
        if temp
            iscue1{ii}(jj) = 1;
        else
            iscue1{ii}(jj) = 2;
        end
    end
end
iscue2 = cell(length(l2), 1);
for ii = 1:length(l2) %classify cue/traj ensembles
    iscue2{ii} = zeros(length(l2{ii}), 1);
    for jj = 1:length(l2{ii})
        if isempty(l2{ii}{jj}); continue; end
        temp = any(s2{ii}{jj} < cue_centres & e2{ii}{jj} > cue_centres & l2{ii}{jj} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
        if temp
            iscue2{ii}(jj) = 1;
        else
            iscue2{ii}(jj) = 2;
        end
    end
end



%% Manual classification
temp = cell(size(iscue2));
figure
for ii = 1:length(iscue2)
    for jj = 1:length(iscue2{ii})
        [~, order] = max(swr_clust_stack{2}{ii}{jj}, [], 1);
        [~, order] = sort(order);
        imagesc(swr_clust_stack{2}{ii}{jj}(:, order)')
        temp{ii}(jj) = input('class: ');
    end
end
iscue2 = temp;


%%
temp = [];
sig = [];
any_sig = cell(length(iscue2), 1);
stack_pairs = {};
for ii = 1:length(iscue2)
    cue = find(iscue2{ii} == 1);
    traj = find(iscue2{ii} == 2);
    any_sig{ii} = false(length(iscue2{ii}), 1);
    for cc = 1:length(cue)
        for tt = 1:length(traj)
            t = ~any(isnan(hiepi_z{2}{ii}), 2);
            t = linspace(0, sum(t) / 19, sum(t));
            [r, ~, p] = bxcorr(hiepi_z{2}{ii}(~any(isnan(hiepi_z{2}{ii}), 2), cue(cc)), hiepi_z{2}{ii}(~any(isnan(hiepi_z{2}{ii}), 2), traj(tt)), t);
            temp = cat(2, temp, r( ceil(length(r)/2)-19 : ceil(length(r)/2)+19 ));
            sig = cat(1, sig, p);
            stack_pairs = cat(1, stack_pairs, [{swr_clust_stack{2}{ii}{cue(cc)}}, {swr_clust_stack{2}{ii}{traj(tt)}}]);
            
            any_sig{ii}(traj(tt)) = any([any_sig{ii}(traj(tt)), any(temp(:, end)>p)]);
        end
    end
end

order = max(temp, [], 1);
[~, order] = sort(order);
figure
imagesc(temp(:, order)');
axis image
colorbar
caxis([0 .3])
colormap jet

r = [];
for ii = 1:size(stack_pairs, 1)
    r = [r; corr(mean(stack_pairs{ii, 1}')', mean(stack_pairs{ii, 2}')')];
end

% idx = find(any(temp > .05, 1));
idx = find(any(temp'>sig, 2));
for ii = 1:length(idx)
    if ~mod(ii - 1, 6)
        figure
    end
    subplot(6, 3, (mod(ii - 1, 6)+1)*3-2);
    [~, order] = max(stack_pairs{idx(ii), 1}, [], 1);
    [~, order] = sort(order);
    imagesc(-stack_pairs{idx(ii), 1}(:, order)');
    colormap bone
    
    subplot(6, 3, (mod(ii - 1, 6)+1)*3-1);
    [~, order] = max(stack_pairs{idx(ii), 2}, [], 1);
    [~, order] = sort(order);
    imagesc(-stack_pairs{idx(ii), 2}(:, order)');
    colormap bone
    
    subplot(6, 3, (mod(ii - 1, 6)+1)*3);
    plot(mean(stack_pairs{idx(ii), 1}'))
    hold on
    plot(mean(stack_pairs{idx(ii), 2}'))
    title(num2str(r(idx(ii))))
end


%%
s2 = [s2{:}];
e2 = [e2{:}];
iscue2 = cell2mat(iscue2);
any_sig = cell2mat(any_sig);

bins = 25;
edges = linspace(0, 150, bins+1);

P_se = accumarray([discretize(cell2mat(s2(iscue2 == 2)'), edges) discretize(cell2mat(e2(iscue2 == 2)'), edges)], 1, [bins bins]);
% P_se = accumarray([discretize(cell2mat(s2((iscue2 == 2) & any_sig)'), edges) discretize(cell2mat(e2((iscue2 == 2) & any_sig)'), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_s_cond_e = P_se ./ sum(P_se, 1);
P_e_cond_s = P_se ./ sum(P_se, 2);

P_e_cond_s = imgaussfilt(P_e_cond_s, 1.5);

figure
imagesc(P_e_cond_s)
rbmap('caxis', [0 1]);
colorbar
caxis([0 1])
axis image


%%
temp_rsc = temp;
save('temp_rsc', 'temp_rsc')
temp_ee = temp;
save('temp_ee', 'temp_ee')


%%
load('/mnt/storage/HaoRan/RRR_motor/M2/temp_rsc.mat')
load('/mnt/storage/HaoRan/RRR_motor/M2/temp_ee.mat')

temp = [temp_rsc, temp_ee];

order = max(temp, [], 1);
[~, order] = sort(order);
figure
imagesc(temp(:, order)');
axis image
colorbar
caxis([0 .3])




