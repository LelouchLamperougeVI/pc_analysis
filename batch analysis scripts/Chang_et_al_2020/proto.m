clear all

root = '/mnt/storage/rrr_magnum/M2';
animal = 'EE003';
date = '2018_12_29';
session = '3';

ass = ensemble(fullfile(root, animal, date, session, [date '_' session '.abf']));

%%
ass.set_ops('e_size',10);
ass.set_ops('clust_method','thres');
ass.cluster;
ass.set_ops('order','pc');

%% plot rest-run
figure
subplot(1, 2, 1);
deconv = ass.analysis.deconv(:, ass.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(-deconv');
colormap gray
subplot(1, 2, 2);
deconv = ass.twop.deconv(:, ass.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(-deconv');
colormap gray


%% Data tree
clear all

date_fcn = @(x) datetime(x, 'inputformat', 'yyyy_MM_dd');
root = '/mnt/storage/rrr_magnum/M2';
% tree = dirTree(root, 'isdir', [1 1 1 1 1 0], 're', {nan, nan, nan, '3', 'plane.', 'deconv.mat'}, 'fcn', {nan,nan,date_fcn,nan,nan,nan});
tree = dirTree(root, 'isdir', [1 1 1 0], 're', {nan, nan, nan, '_2.abf'}, 'fcn', {nan,nan,date_fcn,nan});
anal_list = tree.fullfile;
% anal_list( cellfun(@length, anal_list) < 60 ) = [];

exceptions = {'/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_2.abf';
    '/mnt/storage/rrr_magnum/M2/EE001/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC003/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_20/2018_12_20_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_28/2018_12_28_2.abf'};

anal_list = setdiff(anal_list, exceptions);
clearvars -except anal_list

%% make batch analysis
% for ii = 8:length(anal_list)
for ii = 54:length(anal_list)
    lfp = LFP(anal_list{ii});
    lfp.perform_analysis;
    lfp.save('analysis');
    lfp.save('chan');
end


%% check analysis
clear all

root = '/mnt/storage/rrr_magnum/M2';
tree = dirTree(root, 'isdir', [1 1 1 0], 're', {nan, nan, nan, 'analysis.mat'}, 'fcn', {nan,nan,nan,nan});
anal_list = tree.fullfile;
anal_list( cellfun(@length, anal_list) < 56 ) = [];

for ii = 1:length(anal_list)
    load(anal_list{ii});
    disp([anal_list{ii} ': ' num2str(length(analysis.pc_list) / size(analysis.original_deconv, 2))]);
end


%% Aubrey's animals with mutable cues
stack1 = run1.analysis.stack;
stack2 = run2.analysis.stack;

[~, order] = max(stack1, [], 1);
[~, order] = sort(order);
order = intersect(order, run1.analysis.pc_list, 'stable');

figure
subplot(2, 2, 1);
imagesc(stack1(:, order)')
colormap jet
subplot(2, 2, 2);
imagesc(stack2(:, order)')
colormap jet

[~, order] = max(stack2, [], 1);
[~, order] = sort(order);
order = intersect(order, run2.analysis.pc_list, 'stable');

subplot(2, 2, 3);
imagesc(stack1(:, order)')
colormap jet
subplot(2, 2, 4);
imagesc(stack2(:, order)')
colormap jet


%% Fig1C
cd /home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB

run = LFP('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_2.abf');
run.perform_analysis;

rest1 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_1.abf');
rest2 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');

rest1.import_analysis(run.analysis)
rest1.set_ops('order','pc');
rest2.import_analysis(run.analysis)
rest2.set_ops('order','pc');

%%
figure
ax1(1) = subplot(5, 3, [1 4 7]);
deconv = rest1.twop.deconv(:, rest1.analysis.order);
deconv = fast_smooth(deconv, rest1.ops.sig * rest1.twop.fs);
deconv = (deconv - min(deconv)) ./ range(deconv);
imagesc('xdata', rest1.twop.ts, 'cdata', -deconv');
colormap gray
ax1(2) = subplot(5, 3, 10);
plot(rest1.behavior.ts, rest1.behavior.speed);
ylim([-20 60])
ax1(3) = subplot(5, 3, 13);
plot(rest1.behavior.ts, rest1.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax1, 'x')
xlim([300 400])

ax2(1) = subplot(5, 3, [2 5 8]);
deconv = run.twop.deconv(:, run.analysis.order);
deconv = fast_smooth(deconv, rest1.ops.sig * run.twop.fs);
deconv = (deconv - min(deconv)) ./ range(deconv);
imagesc('xdata', run.twop.ts, 'cdata', -deconv');
colormap gray
ax2(2) = subplot(5, 3, 11);
plot(run.behavior.ts, run.behavior.speed);
ylim([-20 60])
ax2(3) = subplot(5, 3, 14);
plot(run.behavior.ts, run.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax2, 'x')
xlim([200 300])

ax3(1) = subplot(5, 3, [3 6 9]);
deconv = rest2.twop.deconv(:, rest2.analysis.order);
deconv = fast_smooth(deconv, rest2.ops.sig * rest2.twop.fs);
deconv = (deconv - min(deconv)) ./ range(deconv);
imagesc('xdata', rest2.twop.ts, 'cdata', -deconv');
colormap gray
ax3(2) = subplot(5, 3, 12);
plot(rest2.behavior.ts, rest2.behavior.speed);
ylim([-20 60])
ax3(3) = subplot(5, 3, 15);
plot(rest2.behavior.ts, rest2.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax3, 'x')
xlim([300 400])

%% Fig1D
mimg = rest2.topo.mimg;
mimg = (mimg - min(mimg(:))) ./ range(mimg(:));
mimg = mimg ./ .25; %adjust contrast
mimg(mimg>1) = 1;
mask = rest2.topo.maskNeurons;

f = rest2.twop.deconv(rest2.twop.ts > 310.5 & rest2.twop.ts < 311.5, :);
f = mean(f) ./ range(rest2.twop.deconv);
im = repmat(mimg, 1, 1, 3);
c = rbmap('caxis',[0 max(f)]);
f = discretize(f, size(c, 1));
for ii = 1:length(f)
    [x, y] = find(mask == ii);
    for jj = 1:length(x)
        im(x(jj), y(jj), :) = permute(c(f(ii), :), [1 3 2]);
    end
end

imwrite(im, 'fig1d_mask_react.png')
% figure
% imshow(im);

[~, idx] = sort(f, 'descend');
[~, idx] = max( run.analysis.stack(:, idx(1:5)), [], 1 );
idx = round(median(idx));

f = run.analysis.stack(idx-1:idx+1, :);
f = mean(f);
im = repmat(mimg, 1, 1, 3);
c = rbmap('caxis',[0 max(f)]);
f = discretize(f, size(c, 1));
for ii = 1:length(f)
    [x, y] = find(mask == ii);
    for jj = 1:length(x)
        im(x(jj), y(jj), :) = permute(c(f(ii), :), [1 3 2]);
    end
end

imwrite(im, 'fig1d_mask_run.png')
% figure
% imshow(im);

%% Fig2
rest2.set_ops('e_size',5);
rest2.set_ops('clust_method','thres');
rest2.set_ops('sig', .2);
rest2.remove_mvt;
rest2.cluster;
rest2.topography;
rest2.detect_sce;

% rest2.plot('clust_topo')
rest2.set_ops('order','cluster')

rest2.plot('sce')
rest2.plot('clust_corr')
rest2.plot('tree')
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{4})
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{5})
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{7})

plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{4})
plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{5})
plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{7})
% nice clusts: 4, 5, 7

%%
figure
deconv = fast_smooth(rest2.twop.deconv, rest2.ops.sig * rest2.twop.fs * 4);
deconv = (deconv - min(deconv, [], 'omitnan')) ./ range(deconv);
imagesc('xdata', rest2.twop.ts, 'cdata', -deconv(:, rest2.ensembles.clust_order)');
colormap gray
% xlim([355 415])

hold on
[x, y] = find(rest2.twop.deconv(:, rest2.ensembles.clust_order) > std(rest2.twop.deconv(:, rest2.ensembles.clust_order), 'omitnan').*5);

colours = rest2.ensembles.colours;
% colours([4, 6, 1], :) = colours([6, 1, 4], :);

for ii = 1:length(rest2.ensembles.clust)
    idx = any(rest2.ensembles.clust_order(:) == rest2.ensembles.clust{ii}(:)', 2);
    idx = ismember(y, find(idx));
    idx = idx & (rest2.hiepi.z(x, ii) > (std(rest2.hiepi.z(:, ii), 'omitnan') * 5));
    plot(rest2.twop.ts(x(idx)), y(idx), '.', 'MarkerFaceColor', colours(ii, :), 'MarkerSize', 20)
end



%% Batch make/save analysis old M2
chan = [1; 2; 3; 5; NaN; NaN; NaN];

root = '/mnt/storage/HaoRan/RRR_motor/M2';

animals = dir(fullfile(root, 'RSC*'));
animals = {animals.name};

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        disp(['current session: ' fullfile(animals{a}, sessions{s})])
        save(fullfile(root, animals{a}, sessions{s}, 'channels.mat'), 'chan');
        lfp = LFP(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_2.abf']));
        lfp.perform_analysis;
        lfp.save('analysis');
    end
end

%% Batch make/save analysis new M2
clear all

root = '/mnt/storage/rrr_magnum/M2';

animals = dir(fullfile(root, 'E*'));
animals = {animals.name};

for a = 4:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        lfp = LFP(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_2.abf']));
        lfp.perform_analysis;
        lfp.save('analysis');
    end
end


%% Core analysis
clear all

traj_thres = .2;

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
% animals = dir(fullfile(root, 'RSC*'));

root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

EV = [];
frac = []; %[rest1 rest2 run]
frac_clust = cell(2,1);
clust_stacks = cell(2,1);
trajectories = cell(2,1);
loc = {};
loc_clust = cell(2,1);
swr_stack = cell(2,1);
swr_hiepi = cell(2,1);
hiepi_lfp_pw = cell(2,1);
hiepi_z = cell(2,1);
hiepi_psth = cell(2, 1);
masks = cell(2,1);
clusts = cell(2,1);
SI = [];
si_frac = [];
swr_clust_stack = cell(2,1);
whole_stack = {};
pc_list = {};
pca_var = []; % percent explained variance for first 50 components (session X component X rest1/2)
pc_width = {};

file = {}; % list of sessions

ll = cell(2, 1); % log likelihood from pc cc models
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        file = cat(1, file, fullfile(root, animals{a}, sessions{s}));
        
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.set_ops('e_size',5);
        rest1.set_ops('clust_method','thres');
        rest1.set_ops('sig', .2);
        rest1.remove_mvt;
        rest1.cluster;
%         rest1.swr_window;
        rest1.hPICA;
%         rest1.topography;
%         rest1.detect_sce;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
%         rest2.swr_window;
        rest2.hPICA;
%         rest2.topography;
%         rest2.detect_sce;
        [temp1, temp2] = ev(rest1, rest1.analysis.original_deconv, rest2);
        EV = cat(1, EV, [temp1 temp2]);
        
        frac = cat(1, frac, [length(intersect( cell2mat(rest1.ensembles.clust), rest1.analysis.pc_list )) / length(cell2mat(rest1.ensembles.clust)), ...
        length(intersect( cell2mat(rest2.ensembles.clust), rest1.analysis.pc_list )) / length(cell2mat(rest2.ensembles.clust)), ...
        length(rest1.analysis.pc_list) / size(rest1.twop.deconv, 2)]);
        
        si_frac = cat(1, si_frac, {{rest1.analysis.SI(cell2mat(rest1.ensembles.clust))}, {rest1.analysis.SI(cell2mat(rest2.ensembles.clust))}, {rest1.analysis.SI}});
        
        try % RSC animals don't have electrode
            %hPICA
            swr_hiepi{1} = cat(1, swr_hiepi{1}, {rest1.hiepi.swr_react_strength});
            swr_hiepi{2} = cat(1, swr_hiepi{2}, {rest2.hiepi.swr_react_strength});
        catch
        end
        
        hiepi_lfp_pw{1} = cat(1, hiepi_lfp_pw{1}, {rest1.hiepi.reactivations});
        hiepi_lfp_pw{2} = cat(1, hiepi_lfp_pw{2}, {rest2.hiepi.reactivations});
        
        %hPICA z traces
        hiepi_z{1} = cat(1, hiepi_z{1}, {rest1.hiepi.z});
        hiepi_z{2} = cat(1, hiepi_z{2}, {rest2.hiepi.z});
        
        %hPICA z psth
        hiepi_psth{1} = cat(1, hiepi_psth{1}, {rest1.hiepi.psth});
        hiepi_psth{2} = cat(1, hiepi_psth{2}, {rest2.hiepi.psth});
        
        %ROI masks
        masks{1} = cat(1, masks{1}, {rest1.topo.maskNeurons});
        masks{2} = cat(1, masks{2}, {rest2.topo.maskNeurons});
        
        % whole stack for cross days analysis
        whole_stack{end+1} = rest1.analysis.stack;
        pc_list{end+1} = rest1.analysis.pc_list;
        
        % v checkpoint
    
%         swr_stack{1} = cat(1, swr_stack{1}, {rest1.ensembles.swr.all});
%         swr_stack{2} = cat(1, swr_stack{2}, {rest2.ensembles.swr.all});
        clusts{1} = cat(2, clusts{1}, {rest1.ensembles.clust});
        clusts{2} = cat(2, clusts{2}, {rest2.ensembles.clust});
        SI = cat(2, SI, {rest1.analysis.SI});
        
        pc_width = cat(2, pc_width, {rest1.analysis.width});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        swr_clust_stack{1} = cat(2, swr_clust_stack{1}, {cell(1, length(rest1.ensembles.clust))});
        
        % pc_cc characterisation - block 1/3
        analysis = rest2.analysis;
        thres=noRun(analysis.behavior.unit_vel);
        thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
        n = analysis.original_deconv;
        dt = .05; % 300 milliseconds
        dt = round(analysis.fs * dt);
        x = analysis.behavior.unit_pos;
        
        for c = 1:length(rest1.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest1.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            frac_clust{1} = cat(1, frac_clust{1}, [length(list) length(list)/length(rest1.ensembles.clust{c})]);
            clust_stacks{1} = cat(1, clust_stacks{1}, {stack});
            swr_clust_stack{1}{end}{c} = stack;
            trajectories{1} = cat(2, trajectories{1}, any(stack > traj_thres, 2));
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
            
            % pc_cc characterisation - block 2/3
            list = rest1.ensembles.clust{c};
            md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false, 'comp', 'intra');
            ll{1} = cat(1, ll{1}, {[md.pc.l(:), md.cc.l(:)]});
        end
        temp = cell2mat({rest1.analysis.width{~cellfun(@isempty, rest1.analysis.width)}}');
        loc{end+1} =  temp(:, 2);
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{1} = cat(1, loc_clust{1}, {temp});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        swr_clust_stack{2} = cat(2, swr_clust_stack{2}, {cell(1, length(rest2.ensembles.clust))});
        
        for c = 1:length(rest2.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest2.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            frac_clust{2} = cat(1, frac_clust{2}, [length(list) length(list)/length(rest2.ensembles.clust{c})]);
            clust_stacks{2} = cat(1, clust_stacks{2}, {stack});
            swr_clust_stack{2}{end}{c} = stack;
            trajectories{2} = cat(2, trajectories{2}, any(stack > traj_thres, 2));
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end

            % pc_cc characterisation - block 3/3
            list = rest2.ensembles.clust{c};
            md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false, 'comp', 'intra');
            ll{2} = cat(1, ll{2}, {[md.pc.l(:), md.cc.l(:)]});
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{2} = cat(1, loc_clust{2}, {temp});
        
        [~, ~, score1] = pca(rest1.ensembles.R);
        [~, ~, score2] = pca(rest2.ensembles.R);
        
        pca_var = cat(1, pca_var, permute( [score1(1:50) ./ sum(score1), score2(1:50) ./ sum(score2)], [3 1 2] ));
    end
end


%% Fig2e-f
clear all
frac1 = load('/mnt/storage/HaoRan/RRR_motor/M2/frac.mat');
frac_clust1 = load('/mnt/storage/HaoRan/RRR_motor/M2/frac_clust.mat');
frac2 = load('/mnt/storage/rrr_magnum/M2/frac.mat');
frac_clust2 = load('/mnt/storage/rrr_magnum/M2/frac_clust.mat');

frac = cat(1, frac1.frac, frac2.frac);
frac_clust = cell(2, 1);
frac_clust{1} = cat(1, frac_clust1.frac_clust{1}, frac_clust2.frac_clust{1});
frac_clust{2} = cat(1, frac_clust1.frac_clust{2}, frac_clust2.frac_clust{2});

figure
boxplot(frac);
[~,~,stats] = kruskalwallis(frac);
multcompare(stats)

figure
boxplot([frac(:,1) - frac(:,3), frac(:,2) - frac(:,3)]);
ylim([-.5 1])
p = ranksum(frac(:,1) - frac(:,3), frac(:,2) - frac(:,3))
p = signrank(frac(:,1) - frac(:,3), 0, 'tail', 'right')
p = signrank(frac(:,2) - frac(:,3), 0, 'tail', 'right')


%% Fig2g
clear all
si_frac1 = load('/mnt/storage/HaoRan/RRR_motor/M2/si_frac.mat');
si_frac2 = load('/mnt/storage/rrr_magnum/M2/si_frac.mat');
si_frac = cat(1, si_frac1.si_frac, si_frac2.si_frac);

temp = cell(3,1);
for ii = 1:size(si_frac, 1)
    temp{1} = cat(2, temp{1}, si_frac{ii, 1}{1});
    temp{2} = cat(2, temp{2}, si_frac{ii, 2}{1});
    temp{3} = cat(2, temp{3}, si_frac{ii, 3}{1});
end

figure
cdfplot(temp{1})
hold on
cdfplot(temp{2})
cdfplot(temp{3})

figure
boxplot([temp{1} temp{2}], [ones(1, length(temp{1})) 2.*ones(1, length(temp{2}))]);

[~, p] = kstest2(temp{1}, temp{2})
p = ranksum(temp{1}, temp{2})


%% Combine run for EE003/2018_12_19
lfp1 = LFP('/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_2.abf');
behav1 = lfp1.behavior;
lfp2 = LFP('/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_3.abf');
behav2 = lfp2.behavior;
ts1 = lfp1.twop.ts;
ts2 = lfp2.twop.ts;
dec1 = lfp1.twop.deconv;
dec2 = lfp2.twop.deconv;
ts = cat(1, ts1, ts2 + ts1(end) + median(diff(ts1)));
dec = cat(1, dec1, dec2);

beh = behav1;
beh.speed_raw_noSmooth = [beh.speed_raw_noSmooth behav2.speed_raw_noSmooth];
beh.speed_raw = [beh.speed_raw; behav2.speed_raw];

beh.ts = [beh.ts behav2.ts + ts1(end) + median(diff(ts1))];

beh.pos_norm = [beh.pos_norm behav2.pos_norm];
beh.trial= [beh.trial behav2.trial+max(beh.trial)];
beh.speed = [beh.speed behav2.speed];
beh.pos_cum = [beh.pos_cum behav2.pos_cum];
beh.pos_raw = [beh.pos_raw behav2.pos_raw];

tcs.tt = ts';
[beh, dec] = convert_behavior(beh, tcs, dec);
analysis = pc_batch_analysis(beh, dec);


%%
clear all

frac_thres = 3;

frac_clust1 = load('/mnt/storage/HaoRan/RRR_motor/M2/frac_clust.mat');
frac_clust1 = frac_clust1.frac_clust;
frac_clust2 = load('/mnt/storage/rrr_magnum/M2/frac_clust.mat');
frac_clust2 = frac_clust2.frac_clust;

trajectories1 = load('/mnt/storage/HaoRan/RRR_motor/M2/traj.mat');
trajectories1 = trajectories1.trajectories;
trajectories2 = load('/mnt/storage/rrr_magnum/M2/traj.mat');
trajectories2 = trajectories2.trajectories;

trajectories = [trajectories1{1} trajectories2{1}];
frac_clust = [frac_clust1{1}(:,1); frac_clust2{1}(:,1)];
trajectories = trajectories(:, frac_clust >= frac_thres);
[l1, s1, e1] = traj_length(trajectories);
% l = cellfun(@(x) any(x > 20), l);

trajectories = [trajectories1{2} trajectories2{2}];
frac_clust = [frac_clust1{2}(:,1); frac_clust2{2}(:,1)];
trajectories = trajectories(:, frac_clust >= frac_thres);
[l2, s2, e2] = traj_length(trajectories);
% l = cellfun(@(x) any(x > 20), l);

belt;

figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

figure
histogram(cell2mat(s1), 50)
hold on
histogram(cell2mat(e1), 50)
plot(linspace(0, 150, 50), belt_idx.*10)

figure
histogram(cell2mat(s2), 50)
hold on
histogram(cell2mat(e2), 50)
plot(linspace(0, 150, 50), belt_idx.*10)


%% Figure 3
clear all

circ_dist = @(X0, X) min(cat(3, mod(X0 - X, 50), mod(X - X0, 50)), [], 3) .* 150 ./ 50;

thres = 3; %min number of pc per ensemble

loc1 = load('/mnt/storage/HaoRan/RRR_motor/M2/loc.mat');
loc2 = load('/mnt/storage/rrr_magnum/M2/loc.mat');
loc = cat(2, loc1.loc, loc2.loc);
loc_clust{1} = cat(1, loc1.loc_clust{1}, loc2.loc_clust{1});
loc_clust{2} = cat(1, loc1.loc_clust{2}, loc2.loc_clust{2});

d1 = []; d2 = [];
loc_cum1 = []; loc_cum2 = [];

s1 = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = cellfun(@(x) accumarray(x + 1, 1), loc_clust{1}, 'uniformoutput', false);
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 1:length(g{ii})
        d = loc{ii}(loc_clust{1}{ii} == g{ii}(jj));
        loc_cum1 = cat(1, loc_cum1, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d1 = cat(1, d1, mean(d));
    end
end
g = arrayfun(@(x) ismember(loc_clust{1}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s1 = arrayfun(@(x) s1{x}(g{x}), 1:length(s1), 'uniformoutput', false);

s2 = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = cellfun(@(x) accumarray(x + 1, 1), loc_clust{2}, 'uniformoutput', false);
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 1:length(g{ii})
        d = loc{ii}(loc_clust{2}{ii} == g{ii}(jj));
        loc_cum2 = cat(1, loc_cum2, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d2 = cat(1, d2, mean(d));
    end
end
g = arrayfun(@(x) ismember(loc_clust{2}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s2 = arrayfun(@(x) s2{x}(g{x}), 1:length(s2), 'uniformoutput', false);

figure
cdfplot(cell2mat(s1'));
hold on
cdfplot(cell2mat(s2'));

[~, p] = kstest2(cell2mat(s1'), cell2mat(s2'));
disp(['kstest silhouette: ' num2str(p)])

figure
boxplot([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)])

% p = kruskalwallis([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)]);
p = ranksum(cell2mat(s1'), cell2mat(s2'));
disp(['ranksum silhouette: ' num2str(p)])

figure;
violin({d1, d2}, 'labels', {'rest1', 'rest2'}, 'bandwidth', 2)
% violin({d1, d2}, 'labels', {'rest1', 'rest2'}, 'scatter')

% p = kruskalwallis([d1; d2], [ones(length(d1), 1); 2 .* ones(length(d2), 1)]);
p = ranksum(d1, d2);
disp(['ranksum dist: ' num2str(p)])

figure;
histogram(loc_cum1, 50, 'normalization', 'probability')
hold on
histogram(loc_cum2, 50, 'normalization', 'probability')


%% Trajectories
clear all

clust_stacks1 = load('/mnt/storage/HaoRan/RRR_motor/M2/clust_stacks.mat');
clust_stacks2 = load('/mnt/storage/rrr_magnum/M2/clust_stacks.mat');
clust_stacks{1} = cat(1, clust_stacks1.clust_stacks{1}, clust_stacks2.clust_stacks{1});
clust_stacks{2} = cat(1, clust_stacks1.clust_stacks{2}, clust_stacks2.clust_stacks{2});

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
    else
        iscue1(ii) = 2;
    end
end
iscue2 = zeros(length(l2), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l2) %classify cue/traj ensembles
    if isempty(l2{ii}); continue; end
    temp = any(s2{ii} < cue_centres & e2{ii} > cue_centres & l2{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue2(ii) = 1;
    else
        iscue2(ii) = 2;
    end
end

figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

figure
histogram(cell2mat(s1(iscue1 == 2)), 50, 'normalization', 'probability')
title('rest1 start')
xlim([0 150]); ylim([0 .1])
figure
histogram(cell2mat(e1(iscue1 == 2)), 50, 'normalization', 'probability')
title('rest1 end')
xlim([0 150]); ylim([0 .1])
% plot(linspace(0, 150, 50), belt_idx.*.01)

figure
histogram(cell2mat(s2(iscue2 == 2)), 50, 'normalization', 'probability')
title('rest2 start')
xlim([0 150]); ylim([0 .1])
figure
histogram(cell2mat(e2(iscue2 == 2)), 50, 'normalization', 'probability')
title('rest2 end')
xlim([0 150]); ylim([0 .1])
% plot(linspace(0, 150, 50), belt_idx.*.01)


bins = 50;

edges = linspace(0, 150, bins+1);
P_se = accumarray([discretize(cell2mat(s1(iscue1 == 2)), edges) discretize(cell2mat(e1(iscue1 == 2)), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_s_cond_e = P_se ./ sum(P_se, 1);
P_e_cond_s = P_se ./ sum(P_se, 2);

P_e_cond_s(isnan(P_e_cond_s)) = 0;
P_e_cond_s = imgaussfilt(P_e_cond_s, 2);
% figure
% imagesc(P_se)
% colormap jet
% colorbar
% rbmap('caxis', [0 1]);
% axis image
figure
imagesc(P_e_cond_s)
% rbmap('caxis', [0 1]);
colormap jet
colorbar
% caxis([0 .2])
axis image
% figure
% imagesc(P_s_cond_e)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image


P_se = accumarray([discretize(cell2mat(s2(iscue2 == 2)), edges) discretize(cell2mat(e2(iscue2 == 2)), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_s_cond_e = P_se ./ sum(P_se, 1);
P_e_cond_s = P_se ./ sum(P_se, 2);

P_e_cond_s(isnan(P_e_cond_s)) = 0;
P_e_cond_s = imgaussfilt(P_e_cond_s, 2);
% figure
% imagesc(P_se)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image
figure
imagesc(P_e_cond_s)
% rbmap('caxis', [0 1]);
colormap jet
colorbar
% caxis([0 .2])
axis image
% figure
% imagesc(P_s_cond_e)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image


%% Ensemble class pie chart
figure
subplot(1,2,1);
pie(accumarray(iscue1+1, 1), [0 1 1], {'none', 'cue', 'traj'});
subplot(1,2,2);
pie(accumarray(iscue2+1, 1), [0 1 1], {'none', 'cue', 'traj'});

%% Fig3c
figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

for n = [26, 40, 41, 51, 68, 108, 109, 295]
    [~, idx] = max(clust_stacks{2}{n}); 
    [~, idx] = sort(idx);
    figure
    imagesc('xdata', linspace(0, 150, 50), 'cdata', -clust_stacks{2}{n}(:, idx)'); 
    title([num2str(n) '; length:' num2str(l2{n}') '; start:' num2str(s2{n}') '; end:' num2str(e2{n}')]); 
    if iscue2(n) == 1
        xlabel('cue')
    else
        xlabel('traj')
    end
    colormap bone
end
%% All trajectories
traj_idx = find(~cellfun(@isempty, l1));
k = 5;
for n = 1:length(traj_idx)
    if ~mod(n-1, k^2)
        figure
    end
    subplot(k, k, mod(n-1, k^2) + 1);
    [~, idx] = max(clust_stacks{1}{traj_idx(n)}); 
    [~, idx] = sort(idx);
    imagesc('xdata', linspace(0, 150, 50), 'cdata', -clust_stacks{1}{traj_idx(n)}(:, idx)'); 
    title([num2str(traj_idx(n)) '; length:' num2str(l1{traj_idx(n)}') '; start:' num2str(s1{traj_idx(n)}') '; end:' num2str(e1{traj_idx(n)}')]); 
    if iscue1(traj_idx(n)) == 1
        xlabel('cue')
    else
        xlabel('traj')
    end
    colormap bone
end

figure

traj_idx = find(~cellfun(@isempty, l2));
k = 5;
for n = 1:length(traj_idx)
    if ~mod(n-1, k^2)
        figure
    end
    subplot(k, k, mod(n-1, k^2) + 1);
    [~, idx] = max(clust_stacks{2}{traj_idx(n)}); 
    [~, idx] = sort(idx);
    imagesc('xdata', linspace(0, 150, 50), 'cdata', -clust_stacks{2}{traj_idx(n)}(:, idx)'); 
    title([num2str(traj_idx(n)) '; length:' num2str(l2{traj_idx(n)}') '; start:' num2str(s2{traj_idx(n)}') '; end:' num2str(e2{traj_idx(n)}')]); 
    if iscue2(traj_idx(n)) == 1
        xlabel('cue')
    else
        xlabel('traj')
    end
    colormap bone
end

%% Reconstructed trajectories
P_e_cond_s(isnan(P_e_cond_s)) = 0;
P_se(isnan(P_se)) = 0;

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


%% bad LFPs
bad: EC002, EE003
really good: EC007, EE006


%%
clear all
load('/mnt/storage/rrr_magnum/M2/swr_stack.mat');

sig = 1; % smoothing factor

n_stack = []; % ensemble neurons stack
e_stack = []; % ensemble stack
p_stack = []; % population stack

% for s = 1:length(clusts{1})
%     if size(swr_stack{1}{s}, 1) ~= 39; continue; end
%     for e = 1:length(clusts{1}{s})
%         temp = swr_stack{1}{s}(:, clusts{1}{s}{e}, :);
%         temp = mean(temp, 3, 'omitnan');
%         n_stack = cat(2, n_stack, temp);
%         temp = mean(temp, 2);
%         e_stack = cat(2, e_stack, temp);
%     end
%     temp = swr_stack{1}{s};
%     temp = mean(temp, 3, 'omitnan');
%     temp = mean(temp, 2);
%     p_stack = cat(2, p_stack, temp);
% end

for s = 1:length(clusts{2})
    if size(swr_stack{2}{s}, 1) ~= 39; continue; end
    for e = 1:length(clusts{2}{s})
        temp = swr_stack{2}{s}(:, clusts{2}{s}{e}, :);
        temp = mean(temp, 3, 'omitnan');
        n_stack = cat(2, n_stack, temp);
        temp = mean(temp, 2);
        e_stack = cat(2, e_stack, temp);
    end
    temp = swr_stack{2}{s};
    temp = mean(temp, 3, 'omitnan');
    temp = mean(temp, 2);
    p_stack = cat(2, p_stack, temp);
end

% n_stack = n_stack(10:29,:);
% e_stack = e_stack(10:29,:);
% p_stack = p_stack(10:29,:);

[~, idx] = max(e_stack);
[~, idx] = sort(idx);
figure;
temp = e_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('e_stack')
caxis([-.1 .1])

[~, idx] = max(n_stack);
[~, idx] = sort(idx);
figure;
temp = n_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('n_stack')
caxis([-.1 .1])

[~, idx] = max(p_stack);
[~, idx] = sort(idx);
figure;
temp = p_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('p_stack')
caxis([-.1 .1])

% temp = fast_smooth(e_stack, sig);
% h = errorshade(mean(temp,2,'omitnan'), sem(temp, 2));
% 
% temp = fast_smooth(p_stack, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2), 'h', h)
% 
% temp = cellfun(@(x) mean(x, 3, 'omitnan'), swr_stack{1}, 'uniformoutput', false);
% temp = cell2mat(temp(setdiff(1:68, 35))');
% temp = fast_smooth(temp, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2))
% 
% temp = cellfun(@(x) mean(x, 3, 'omitnan'), swr_stack{2}, 'uniformoutput', false);
% temp = cell2mat(temp(setdiff(1:68, 35))');
% temp = fast_smooth(temp, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2))


%% Fig5a
clear all

animal = 'EC007';
date = '2019_05_29';

lfp = ensemble(fullfile('/mnt/storage/rrr_magnum/M2/', animal, date, [date '_3.abf']));
lfp.set_ops('e_size',5);
lfp.set_ops('clust_method','thres');
lfp.set_ops('sig', .2);
lfp.remove_mvt;
lfp.cluster;
lfp.set_ops('order','cluster')
lfp.detect_sce;

% lfp.swr_window;
% lfp.detect_sce;
% lfp.sce_spectrum;

deconv = lfp.twop.deconv(:, lfp.ensembles.clust_order);
deconv = fast_smooth(deconv, lfp.ops.sig * lfp.twop.fs);
deconv = (deconv - nanmin(deconv)) ./ range(deconv);

figure
ax(1) = subplot(7,1,1:4);
imagesc('xdata', lfp.twop.ts, 'cdata', -deconv')
colormap gray
ax(2) = subplot(7,1,5);
plot(lfp.twop.ts, lfp.ensembles.MUA);
ax(3) = subplot(7,1,6);
plot(lfp.lfp.ts, lfp.lfp.lfp);
hold on
plot(lfp.lfp.swr.swr_on, max(lfp.lfp.lfp) .* ones(length(lfp.lfp.swr.swr_on), 1), 'k*');
ax(4) = subplot(7,1,7);
plot(lfp.lfp.ts, lfp.lfp.swr.swr);
linkaxes(ax, 'x')
xlim([560 620])

deconv = lfp.twop.deconv(:, lfp.ensembles.clust_order);
deconv = fast_smooth(deconv, 0.05 * lfp.twop.fs);
deconv = (deconv - nanmin(deconv)) ./ range(deconv);

figure
ax(1) = subplot(7,1,1:4);
imagesc('xdata', lfp.twop.ts, 'cdata', -deconv')
colormap gray
ax(2) = subplot(7,1,5);
plot(lfp.twop.ts, lfp.ensembles.MUA);
ax(3) = subplot(7,1,6);
plot(lfp.lfp.ts, lfp.lfp.lfp);
hold on
plot(lfp.lfp.swr.swr_on, max(lfp.lfp.lfp) .* ones(length(lfp.lfp.swr.swr_on), 1), 'k*');
ax(4) = subplot(7,1,7);
plot(lfp.lfp.ts, lfp.lfp.swr.swr);
linkaxes(ax, 'x')
xlim([573 575])

figure
ax(1) = subplot(7,1,1:4);
imagesc('xdata', lfp.twop.ts, 'cdata', -deconv')
colormap gray
ax(2) = subplot(7,1,5);
plot(lfp.twop.ts, lfp.ensembles.MUA);
ax(3) = subplot(7,1,6);
plot(lfp.lfp.ts, lfp.lfp.lfp);
hold on
plot(lfp.lfp.swr.swr_on, max(lfp.lfp.lfp) .* ones(length(lfp.lfp.swr.swr_on), 1), 'k*');
ax(4) = subplot(7,1,7);
plot(lfp.lfp.ts, lfp.lfp.swr.swr);
linkaxes(ax, 'x')
xlim([602 604])

%% Fig5c - SCE onset events triggered average CWT spectrogram
% sce = cat(1, lfp.ensembles.clust_SCE(1).SCE.on, lfp.ensembles.clust_SCE(2).SCE.on, lfp.ensembles.clust_SCE(2).SCE.on);
sce = lfp.ensembles.SCE.on;
eta_cwt(lfp.lfp.lfp, lfp.lfp.fs, .2, sce, 'nans', isnan(lfp.lfp.swr.swr_env), 'plot', true, 'nvc', 20);


%% Fig5b - swr onset events triggered average CWT spectrogram
clear all
animal = 'EC007';
date = '2019_05_29';
lfp = LFP(fullfile('/mnt/storage/rrr_magnum/M2/', animal, date, [date '_3.abf']));

[s, t, f] = eta_cwt(lfp.lfp.lfp, lfp.lfp.fs, .2, lfp.lfp.swr.swr_on, 'nans', isnan(lfp.lfp.swr.swr_env), 'plot', true, 'nvc', 20);
ylim([0 300])


%% SWR trajectories
clear all
% load('/mnt/storage/rrr_magnum/M2/swr_stack_2s.mat');
% load('/mnt/storage/rrr_magnum/M2/hiepi.mat');
load('/mnt/storage/rrr_magnum/M2/hiepi3.mat');

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


%% SWR cue/traj
sigma = 1;

% rest 1
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l1)
    for jj = 1:length(l1{ii})
        temp = mean( mean( swr_stack{1}{ii}(:, clusts{1}{ii}{jj}, :), 3 ), 2 );
        if iscue1{ii}(jj) == 2
            traj_ens = cat(2, traj_ens, temp);
        elseif iscue1{ii}(jj) == 1
            cue_ens = cat(2, cue_ens, temp);
        else
            try
                none_ens = cat(2, none_ens, temp);
            catch
            end
        end
    end
end

t = linspace(-2, 2, size(traj_ens, 1));
% temp = fast_smooth(traj_ens, sigma);
temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
% temp = fast_smooth(cue_ens, sigma);
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
% temp = fast_smooth(none_ens, sigma);
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])


% rest 2
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l2)
    for jj = 1:length(l2{ii})
        temp = mean( mean( swr_stack{2}{ii}(:, clusts{2}{ii}{jj}, :), 3 ), 2 );
        if iscue2{ii}(jj) == 2
            traj_ens = cat(2, traj_ens, temp);
        elseif iscue2{ii}(jj) == 1
            cue_ens = cat(2, cue_ens, temp);
        else
            try
                none_ens = cat(2, none_ens, temp);
            catch
            end
        end
    end
end

temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])


%% SWR cue/traj with hiepi
sigma = 1;

% rest 1
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l1)
    for jj = 1:length(l1{ii})
        temp = mean(swr_hiepi{1}{ii}{jj}, 'omitnan')';
        if iscue1{ii}(jj) == 2
            traj_ens = cat(2, traj_ens, temp);
        elseif iscue1{ii}(jj) == 1
            cue_ens = cat(2, cue_ens, temp);
        else
            try
                none_ens = cat(2, none_ens, temp);
            catch
            end
        end
    end
end

t = linspace(-4, 4, size(traj_ens, 1));
idx = t > -1 & t < 1;
% temp = fast_smooth(traj_ens, sigma);
temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
% temp = fast_smooth(cue_ens, sigma);
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
% temp = fast_smooth(none_ens, sigma);
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])

figure
[r, delay, p] = bxcorr(mean(traj_ens(idx, :), 2, 'omitnan'), mean(cue_ens(idx, :), 2, 'omitnan'), t);
f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
plot(f, delay, r, 'k')
yline(p)
text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
title('traj x cue')
xlim([-1 1])
ylim([-1.2 1.2])
legend('xcorr', 'fit')

% figure
% ax(1) = subplot(1, 2, 1);
% [r, delay, p] = bxcorr(mean(none_ens(idx, :), 2, 'omitnan'), mean(traj_ens(idx, :), 2, 'omitnan'), t);
% f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
% plot(f, delay, r, 'k')
% yline(p)
% text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
% title('none x traj')
% ax(2) = subplot(1, 2, 2);
% [r, delay, p] = bxcorr(mean(none_ens(idx, :), 2, 'omitnan'), mean(cue_ens(idx, :), 2, 'omitnan'), t);
% f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
% plot(f, delay, r, 'k')
% yline(p)
% text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
% title('none x cue')
% linkaxes(ax, 'xy');
% xlim([-1 1])
% legend('xcorr', 'fit')


% rest 2
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l2)
    for jj = 1:length(l2{ii})
        temp = mean(swr_hiepi{2}{ii}{jj}, 'omitnan')';
        if iscue2{ii}(jj) == 2
            traj_ens = cat(2, traj_ens, temp);
        elseif iscue2{ii}(jj) == 1
            cue_ens = cat(2, cue_ens, temp);
        else
            try
                none_ens = cat(2, none_ens, temp);
            catch
            end
        end
    end
end

temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])

figure
[r, delay, p] = bxcorr(mean(traj_ens(idx, :), 2, 'omitnan'), mean(cue_ens(idx, :), 2, 'omitnan'), t);
f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
plot(f, delay, r, 'k')
yline(p)
text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
title('traj x cue')
xlim([-1 1])
ylim([-1.2 1.2])
legend('xcorr', 'fit')


%% LFP power cue/traj with hiepi2
sigma = 100;

% rest 1
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l1)
    for jj = 1:length(l1{ii})
        temp = mean(hiepi_lfp_pw{1}{ii}{jj}.swr_pw, 'omitnan')';
%         temp = mean(hiepi_lfp_pw{1}{ii}{jj}.delta_pw, 'omitnan')';
        if iscue1{ii}(jj) == 2
            traj_ens = cat(2, traj_ens, temp);
        elseif iscue1{ii}(jj) == 1
            cue_ens = cat(2, cue_ens, temp);
        else
            try
                none_ens = cat(2, none_ens, temp);
            catch
            end
        end
    end
end

t = linspace(-4, 4, size(traj_ens, 1));
idx = t > -1 & t < 1;
% temp = fast_smooth(traj_ens, sigma);
temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
% temp = fast_smooth(cue_ens, sigma);
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
% temp = fast_smooth(none_ens, sigma);
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])

figure
[r, delay, p] = bxcorr(mean(traj_ens(idx, :), 2, 'omitnan'), mean(cue_ens(idx, :), 2, 'omitnan'), t);
f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
plot(f, delay, r, 'k')
yline(p)
text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
title('traj x cue')
xlim([-1 1])
ylim([-1.2 1.2])
legend('xcorr', 'fit')


% rest 2
cue_ens = [];
traj_ens = [];
none_ens = [];

for ii = 1:length(l2)
    for jj = 1:length(l2{ii})
        temp = mean(hiepi_lfp_pw{2}{ii}{jj}.swr_pw, 'omitnan')';
%         temp = mean(hiepi_lfp_pw{2}{ii}{jj}.delta_pw, 'omitnan')';
        try
            if iscue2{ii}(jj) == 2
                traj_ens = cat(2, traj_ens, temp);
            elseif iscue2{ii}(jj) == 1
                cue_ens = cat(2, cue_ens, temp);
            else
                try
                    none_ens = cat(2, none_ens, temp);
                catch
                end
            end
        catch
            continue
        end
    end
end

temp = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
h = errorshade(t, mean(temp, 2), sem(temp, 2), 'colour', 'b');
hold on
temp = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'r');
temp = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
errorshade(t, mean(temp, 2), sem(temp, 2), 'h', h, 'colour', 'k');
xline(0);
yline(0);
% ylim([-.05 .05])

figure
[r, delay, p] = bxcorr(mean(traj_ens(idx, :), 2, 'omitnan'), mean(cue_ens(idx, :), 2, 'omitnan'), t);
f = fit(delay(delay > -1 & delay < 1), r(delay > -1 & delay < 1), 'gauss1');
plot(f, delay, r, 'k')
yline(p)
text(f.b1, f(f.b1) + .02, [num2str(f.b1) newline '\downarrow'], 'HorizontalAlignment', 'center');
title('traj x cue')
xlim([-1 1])
ylim([-1.2 1.2])
legend('xcorr', 'fit')


%% Create data table for R
temp1 = fast_smooth((none_ens - mean(none_ens)) ./ std(none_ens), sigma);
temp2 = fast_smooth((traj_ens - mean(traj_ens)) ./ std(traj_ens), sigma);
temp3 = fast_smooth((cue_ens - mean(cue_ens)) ./ std(cue_ens), sigma);

t_partitions = [-inf -1;   %before
                -.3  0;  %early
                0    .4    %during
                1    inf]; %after
% t_partitions = [-inf -1;   %before
%                 -.4  .2;  %early
%                 1    inf]; %after
t_partitions = cat(2, t_partitions(:,1) < permute(t, [1 3 2]), t_partitions(:,2) > permute(t, [1 3 2]));
t_partitions = squeeze(and(t_partitions(:,1,:), t_partitions(:,2,:)));

responses = [arrayfun(@(x) temp1(t_partitions(x, :), :), 1:size(t_partitions,1), 'uniformoutput', false);
            arrayfun(@(x) temp2(t_partitions(x, :), :), 1:size(t_partitions,1), 'uniformoutput', false);
            arrayfun(@(x) temp3(t_partitions(x, :), :), 1:size(t_partitions,1), 'uniformoutput', false)];
responses = cellfun(@mean, responses, 'uniformoutput', false);
        
mu_res = cellfun(@(x) mean(x, 'all'), responses);
err_res = cellfun(@(x) sem(x(:)), responses);

figure
bar(mu_res(:));
hold on
errorbar(1:numel(responses), mu_res(:), err_res(:), 'linestyle', 'none');

[gt, gg] = meshgrid(1:size(responses,2), 1:size(responses,1));
gg = arrayfun(@(x, y) y .* ones(size(x{1})), responses, gg, 'uniformoutput', false);
gt = arrayfun(@(x, y) y .* ones(size(x{1})), responses, gt, 'uniformoutput', false);

temp = cellfun(@(x) size(x, 2), responses);
temp = [zeros(1, size(temp,2)); cumsum(temp(1:end-1, :))] + 1;
subjects = arrayfun(@(x, y) meshgrid(y:y+size(x{1},2)-1, 1:size(x{1},1)), responses, temp, 'uniformoutput', false);

y = []; GG = []; GT = []; S = [];
for ii = 1:numel(responses)
    y = cat(1, y, responses{ii}(:));
    GG = cat(1, GG, gg{ii}(:));
    GT = cat(1, GT, gt{ii}(:));
    S = cat(1, S, subjects{ii}(:));
end

% figure
% [~, ~, stats] = anovan(y, {S, GG, GT}, 'random', 1, 'varnames', {'subject', 'ens type', 'time'}, 'model', 'full');
% [c, m, h, gnames] = multcompare(stats, 'dimension', [2 3]);

tbl = table(y, categorical(S), categorical(GG), categorical(GT), 'variablenames', {'response', 'subject', 'type', 'time'});

writetable(tbl, 'test.csv');




%% New data table with neurons, session

% rest 1
stack = [];
type = [];
neur_id = [0];
ens_id = [0];
sess_id = [];

for ii = 1:length(l1)
    for jj = 1:length(l1{ii})
        temp = mean( swr_stack{1}{ii}(:, clusts{1}{ii}{jj}, :), 3 );
        temp = (temp - mean(temp)) ./ std(temp);
        if 77 ~= size(temp,1)
            continue
        end
        stack = cat(2, stack, temp);
        type = cat(1, type, repmat(iscue1{ii}(jj), [size(temp,2) 1]));
        neur_id = cat(1, neur_id, neur_id(end) + (1:size(temp,2))');
        ens_id = cat(1, ens_id, ens_id(end) + ones([size(temp,2) 1]));
        sess_id = cat(1, sess_id, repmat(ii, [size(temp,2) 1]));
    end
end
neur_id(1) = [];
ens_id(1) = [];

t = linspace(-2, 2, size(stack, 1));
t_partitions = [-inf -1;   %before
                -.3  0;  %early
                0    .4    %during
                1    inf]; %after
t_partitions = cat(2, t_partitions(:,1) < permute(t, [1 3 2]), t_partitions(:,2) > permute(t, [1 3 2]));
t_partitions = squeeze(and(t_partitions(:,1,:), t_partitions(:,2,:)));

responses = arrayfun(@(x) stack(t_partitions(x, :), :), 1:size(t_partitions,1), 'uniformoutput', false);
responses = cellfun(@mean, responses, 'uniformoutput', false);
responses = cell2mat(responses')';

figure
temp1 = [mean(responses(type == 0, :)); mean(responses(type == 2, :)); mean(responses(type == 1, :))];
bar(temp1(:));
hold on
temp2 = [sem(responses(type == 0, :)); sem(responses(type == 2, :)); sem(responses(type == 1, :))];
errorbar(1:numel(temp1), temp1(:), temp2(:))

tbl = table(categorical(neur_id), categorical(ens_id), categorical(sess_id), categorical(type), responses(:,1), responses(:,2), responses(:,3), responses(:,4),...
    'variablenames', {'neuron', 'ensemble', 'session', 'type', 't0', 't1', 't2', 't3'});

writetable(tbl, '/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB/r1.csv');

% rest 2
stack = [];
type = [];
neur_id = [0];
ens_id = [0];
sess_id = [];

for ii = 1:length(l2)
    for jj = 1:length(l2{ii})
        temp = mean( swr_stack{2}{ii}(:, clusts{2}{ii}{jj}, :), 3 );
        temp = (temp - mean(temp)) ./ std(temp);
        if 77 ~= size(temp,1)
            continue
        end
        stack = cat(2, stack, temp);
        type = cat(1, type, repmat(iscue2{ii}(jj), [size(temp,2) 1]));
        neur_id = cat(1, neur_id, neur_id(end) + (1:size(temp,2))');
        ens_id = cat(1, ens_id, ens_id(end) + ones([size(temp,2) 1]));
        sess_id = cat(1, sess_id, repmat(ii, [size(temp,2) 1]));
    end
end
neur_id(1) = [];
ens_id(1) = [];

responses = arrayfun(@(x) stack(t_partitions(x, :), :), 1:size(t_partitions,1), 'uniformoutput', false);
responses = cellfun(@mean, responses, 'uniformoutput', false);
responses = cell2mat(responses')';

figure
temp1 = [mean(responses(type == 0, :)); mean(responses(type == 2, :)); mean(responses(type == 1, :))];
bar(temp1(:));
hold on
temp2 = [sem(responses(type == 0, :)); sem(responses(type == 2, :)); sem(responses(type == 1, :))];
errorbar(1:numel(temp1), temp1(:), temp2(:))

tbl = table(categorical(neur_id), categorical(ens_id), categorical(sess_id), categorical(type), responses(:,1), responses(:,2), responses(:,3), responses(:,4),...
    'variablenames', {'neuron', 'ensemble', 'session', 'type', 't0', 't1', 't2', 't3'});

writetable(tbl, '/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB/r2.csv');


%% Supplementary figures

%% supplementary fig1 - lesion vs control from Esteves et al. 2021 jneuro
clear all

% % make them damn objects
% root = '/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation';
% animal = dir(root);
% 
% abfs = {};
% for a = 3:length(animal)
%     date = dir(fullfile(root, animal(a).name));
%     for d = 3:length(date)
%         session = dir(fullfile(root, animal(a).name, date(d).name, '*.abf'));
%         abfs = cat(2, abfs, fullfile(root, animal(a).name, date(d).name, {session.name}));
%     end
% end
% abfs = abfs';
% 
% for ii = 1:length(abfs)
%     obj(ii) = LFP(abfs{ii}, 'planes', 1:3);
%     obj(ii).rm_redund;
%     obj(ii).perform_analysis;
% end


clear all
load('/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation/obj.mat');

% partition analysis structs
l_analysis = [];
c_analysis = [];
for ii = 1:length(obj)
    temp = regexpi(obj(ii).session.animal, '^hl*');
    if isempty(temp)
        c_analysis = [c_analysis; obj(ii).analysis];
    else
        l_analysis = [l_analysis; obj(ii).analysis];
    end
end

% panel a - stacks
c_stack = arrayfun(@(x) x.stack(:, x.pc_list), c_analysis, 'UniformOutput', false);
l_stack = arrayfun(@(x) x.stack(:, x.pc_list), l_analysis, 'UniformOutput', false);
c_stack = cat(2, c_stack{:});
l_stack = cat(2, l_stack{:});

figure
subplot(1, 2, 1);
[~, idx] = max(c_stack, [], 1);
[~, idx] = sort(idx);
imagesc(-c_stack(:, idx)')
rbmap('caxis', [-1 0], 'interp', 255);

subplot(1, 2, 2);
[~, idx] = max(l_stack, [], 1);
[~, idx] = sort(idx);
imagesc(l_stack(:, idx)')
rbmap('caxis', [0 1], 'interp', 255);

% panel b - place field location densities

c_loc = arrayfun(@(x) cat(1, x.width{x.pc_list}), c_analysis, 'UniformOutput', false);
c_loc = cat(1, c_loc{:});
l_loc = arrayfun(@(x) cat(1, x.width{x.pc_list}), l_analysis, 'UniformOutput', false);
l_loc = cat(1, l_loc{:});

g = repelem((1:2)', [size(c_loc, 1); size(l_loc, 1)]);
loc = cat(1, c_loc(:, 2), l_loc(:, 2));

g = gramm('x', loc, 'color', g);
g.stat_bin('nbins', 50, 'geom', 'line', 'fill', 'transparent', 'normalization', 'probability');
g.axe_property('ylim', [0 .08]);
figure
g.draw();

% panel c - number of place fields

c_num_pf = arrayfun(@(x) cellfun(@(y) size(y, 1), x.width(x.pc_list)), c_analysis, 'UniformOutput', false);
c_num_pf = cat(2, c_num_pf{:});
l_num_pf = arrayfun(@(x) cellfun(@(y) size(y, 1), x.width(x.pc_list)), l_analysis, 'UniformOutput', false);
l_num_pf = cat(2, l_num_pf{:});

g = repelem((1:2)', [length(c_num_pf); length(l_num_pf)]);
num_pf = cat(2, c_num_pf, l_num_pf);

figure
g = gramm('x', num_pf, 'color', g);
g.stat_bin('normalization', 'probability', 'edges', .5:1:3.5);
g.axe_property('ylim', [0 1]);
figure
g.draw();

p = ranksum(l_num_pf, c_num_pf)



%% supplementary fig3 - characterisation of resting state dynamics
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat');

% panel a - num neurons per ensemble
frac1 = cat(2, arrayfun(@(x) cellfun(@length, rsc.clusts{1}{x}) ./ size(rsc.whole_stack{x}, 2), 1:length(rsc.whole_stack), 'UniformOutput', false), ...
            arrayfun(@(x) cellfun(@length, ee.clusts{1}{x}) ./ size(ee.whole_stack{x}, 2), 1:length(ee.whole_stack), 'UniformOutput', false));
frac1 = [frac1{:}];
frac2 = cat(2, arrayfun(@(x) cellfun(@length, rsc.clusts{2}{x}) ./ size(rsc.whole_stack{x}, 2), 1:length(rsc.whole_stack), 'UniformOutput', false), ...
            arrayfun(@(x) cellfun(@length, ee.clusts{2}{x}) ./ size(ee.whole_stack{x}, 2), 1:length(ee.whole_stack), 'UniformOutput', false));
frac2 = [frac2{:}];

figure
cdfplot(frac1)
hold on
cdfplot(frac2)
[~, pval] = kstest2(frac1, frac2)
ranksum(frac1, frac2)

% panel b - scatter/corr #ens vs #neurons
num_neur = cat(1, cellfun(@(x) size(x, 2), rsc.whole_stack)', cellfun(@(x) size(x, 2), ee.whole_stack)');
num_ens = [cat(1, cellfun(@length, rsc.clusts{1})', cellfun(@length, ee.clusts{1})'), cat(1, cellfun(@length, rsc.clusts{2})', cellfun(@length, ee.clusts{2})')];
figure
scatter(num_neur, num_ens(:, 2), '.');

hold on
md = fitlm(num_neur, num_ens(:, 1));
plot(num_neur, [ones(length(num_neur), 1), num_neur] * md.Coefficients.Estimate);

% xlim([100 600]);
axis square
title(['r sqr = ' num2str(md.Rsquared.Ordinary) '; r sqr adj = ' num2str(md.Rsquared.Adjusted) '; pval = ' num2str(md.Coefficients.pValue(end))])

% panel c - norm. ensembles per session
figure
boxplot(num_ens ./ num_neur);
ylim([0 .05])
signrank(diff(num_ens, 1, 2) ./ num_neur, 0, 'tail', 'right')

% panel d - rate of reactivation
fs = ee.rest1.twop.fs;

rr1 = [];
for ii = 1:length(rsc.hiepi_lfp_pw{1})
    for jj = 1:length(rsc.hiepi_lfp_pw{1}{ii})
        rr1 = [rr1; length(rsc.hiepi_lfp_pw{1}{ii}{jj}.on) ./ sum(~isnan(rsc.hiepi_z{1}{ii}(:, jj))) .* fs];
    end
end
for ii = 1:length(ee.hiepi_lfp_pw{1})
    for jj = 1:length(ee.hiepi_lfp_pw{1}{ii})
        rr1 = [rr1; length(ee.hiepi_lfp_pw{1}{ii}{jj}.on) ./ sum(~isnan(ee.hiepi_z{1}{ii}(:, jj))) .* fs];
    end
end

rr2 = [];
sparsity2 = [];
for ii = 1:length(rsc.hiepi_lfp_pw{2})
    for jj = 1:length(rsc.hiepi_lfp_pw{2}{ii})
        rr2 = [rr2; length(rsc.hiepi_lfp_pw{2}{ii}{jj}.on) ./ sum(~isnan(rsc.hiepi_z{2}{ii}(:, jj))) .* fs];
    end
end
for ii = 1:length(ee.hiepi_lfp_pw{2})
    for jj = 1:length(ee.hiepi_lfp_pw{2}{ii})
        rr2 = [rr2; length(ee.hiepi_lfp_pw{2}{ii}{jj}.on) ./ sum(~isnan(ee.hiepi_z{2}{ii}(:, jj))) .* fs];
    end
end

pval = ranksum(rr1, rr2)
[~, pval] = kstest2(rr1, rr2)

g = gramm('x', [rr1; rr2], 'color', repelem([1; 2], [length(rr1); length(rr2)]));
g.stat_bin('geom', 'line', 'normalization', 'probability', 'nbins', 25, 'fill', 'transparent');
g.geom_vline('xintercept', [median(rr1), median(rr2)], 'style', '-');

g.axe_property('xlim', [0 .4], 'ylim', [0 .15]);
g.set_title(['Mann-Whitney two-tailed U-test p = ' num2str(pval)]);

figure
g.draw;

figure
cdfplot(rr1)
hold on
cdfplot(rr2)

% panel e - pca explained variance
t = cat(1, rsc.pca_var, ee.pca_var);
% t = diff(t, 1, 3);
t = t(:, 1:5, :);
x = repelem(1:5, size(t, 1), 1, 2);
g = repelem(permute((1:2), [1 3 2]), size(t, 1), size(t, 2), 1);

g = gramm('x', x, 'y', t(:), 'color', g(:));
g.stat_summary('type', 'sem');
g.axe_property('xlim', [1 5], 'ylim', [0 .16]);

t = cat(1, rsc.pca_var, ee.pca_var);
t = diff(t, 1, 3);
t = t(:, 1:5, :);
x = repelem(1:5, size(t, 1), 1);
g(1, 2) = gramm('x', x, 'y', t(:));
g(1, 2).stat_summary('type', 'sem');
g(1, 2).axe_property('xlim', [1 5], 'ylim', [-.03 .03]);
figure
g.draw()


t = cat(1, rsc.pca_var, ee.pca_var);
t = diff(t, 1, 3);

figure
for ii = 1:10
    subplot(4, 4, ii);
    histogram(t(:, ii), 25);
end

t = array2table(t);
t = [array2table((1:size(t, 1))', 'variablename', {'session'}), t];

md = fitrm(t, 't1-t5 ~ 1');

md.ranova
md.mauchly
md.epsilon

md.multcompare('Time', 'comparisontype', 'bonferroni')



%% supp fig 4 - Hypergeometric test (aka Fisher's exact test) for fraction of place cells
clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat');

N = sum(cellfun(@(x) size(x, 2), rsc.whole_stack)) + sum(cellfun(@(x) size(x, 2), ee.whole_stack));
K = sum(cellfun(@length, rsc.pc_list)) + sum(cellfun(@length, ee.pc_list));
k = sum(rsc.frac_clust{1}(:, 1)) + sum(ee.frac_clust{1}(:, 1));
n = sum(rsc.frac_clust{1}(:, 1) ./ rsc.frac_clust{1}(:, 2), 'omitnan') + sum(ee.frac_clust{1}(:, 1) ./ ee.frac_clust{1}(:, 2), 'omitnan');

m = [k,     n - k
     K - k, N - K + k - n];
[~, p] = fishertest(m, 'tail', 'right')

N = sum(cellfun(@(x) size(x, 2), rsc.whole_stack)) + sum(cellfun(@(x) size(x, 2), ee.whole_stack));
K = sum(cellfun(@length, rsc.pc_list)) + sum(cellfun(@length, ee.pc_list));
k = sum(rsc.frac_clust{2}(:, 1)) + sum(ee.frac_clust{2}(:, 1));
n = sum(rsc.frac_clust{2}(:, 1) ./ rsc.frac_clust{2}(:, 2), 'omitnan') + sum(ee.frac_clust{2}(:, 1) ./ ee.frac_clust{2}(:, 2), 'omitnan');

m = [k,     n - k
     K - k, N - K + k - n];
[~, p] = fishertest(m, 'tail', 'right')


% individual ensembles
rest = 1;
whole_stack = cat(2, rsc.whole_stack, ee.whole_stack);
pc_list = cat(2, rsc.pc_list, ee.pc_list);
clusts{1} = cat(2, rsc.clusts{1}, ee.clusts{1});
clusts{2} = cat(2, rsc.clusts{2}, ee.clusts{2});

N = repelem(cellfun(@(x) size(x, 2), whole_stack), cellfun(@length, clusts{rest}));
K = repelem(cellfun(@length, pc_list), cellfun(@length, clusts{rest}));
k = cat(1, rsc.frac_clust{rest}(:, 1), ee.frac_clust{rest}(:, 1));
n = k ./ cat(1, rsc.frac_clust{rest}(:, 2), ee.frac_clust{rest}(:, 2));

N(k < 1) = [];
K(k < 1) = [];
n(k < 1) = [];
k(k < 1) = [];

N = permute(N(:), [2 3 1]);
K = permute(K(:), [2 3 1]);
n = permute(n(:), [2 3 1]);
k = permute(k(:), [2 3 1]);

m = [k,     n - k
     K - k, N - K + k - n];
m = permute(round(m), [2 1 3]);

[~, p1] = arrayfun(@(x) fishertest(m(:, :, x), 'tail', 'right'), 1:size(m, 3));

rest = 2;
whole_stack = cat(2, rsc.whole_stack, ee.whole_stack);
pc_list = cat(2, rsc.pc_list, ee.pc_list);
clusts{1} = cat(2, rsc.clusts{1}, ee.clusts{1});
clusts{2} = cat(2, rsc.clusts{2}, ee.clusts{2});

N = repelem(cellfun(@(x) size(x, 2), whole_stack), cellfun(@length, clusts{rest}));
K = repelem(cellfun(@length, pc_list), cellfun(@length, clusts{rest}));
k = cat(1, rsc.frac_clust{rest}(:, 1), ee.frac_clust{rest}(:, 1));
n = k ./ cat(1, rsc.frac_clust{rest}(:, 2), ee.frac_clust{rest}(:, 2));

N(k < 1) = [];
K(k < 1) = [];
n(k < 1) = [];
k(k < 1) = [];

N = permute(N(:), [2 3 1]);
K = permute(K(:), [2 3 1]);
n = permute(n(:), [2 3 1]);
k = permute(k(:), [2 3 1]);

m = [k,     n - k
     K - k, N - K + k - n];
m = permute(round(m), [2 1 3]);

[~, p2] = arrayfun(@(x) fishertest(m(:, :, x), 'tail', 'right'), 1:size(m, 3));

figure
cdfplot(log(p1))
hold on
cdfplot(log(p2))
xlim([-15 0])

figure
boxplot([log(p1), log(p2)], repelem(1:2, [length(p1), length(p2)]), 'plotstyle', 'compact');
ylim([-15 0])

[~, p] = kstest2(log(p1(log(p1)<-3)), log(p2(log(p2)<-3)), 'tail', 'smaller')
ranksum(log(p1(log(p1)<-3)), log(p2(log(p2)<-3)), 'tail', 'right')


%% supp. fig. 5 - cue/traj PC characterisation
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
    else
        iscue1(ii) = 2;
    end
end
iscue2 = zeros(length(l2), 1); % 1:iscue 2:istraj 0:aint shit
for ii = 1:length(l2) %classify cue/traj ensembles
    if isempty(l2{ii}); continue; end
    temp = any(s2{ii} < cue_centres & e2{ii} > cue_centres & l2{ii} < l_thres, 'all'); % cue centre within traj and length smaller than l_thres
    if temp
        iscue2(ii) = 1;
    else
        iscue2(ii) = 2;
    end
end

clearvars -except rsc ee iscue1 iscue2

width = cat(2, rsc.pc_width, ee.pc_width);
pc_list = cat(2, rsc.pc_list, ee.pc_list);
clusts1 = cat(2, rsc.clusts{1}, ee.clusts{1});
clusts2 = cat(2, rsc.clusts{2}, ee.clusts{2});

num_pf = [];
loc = {};
wid = [];
frac = [];
for ii = 1:length(clusts2)
    for jj = 1:length(clusts2{ii})
        frac = cat(1, frac, length(intersect(clusts2{ii}{jj}, pc_list{ii})) ./ length(clusts2{ii}{jj}));
        temp = width{ii}(intersect(clusts2{ii}{jj}, pc_list{ii}));
        num_pf = cat(1, num_pf, mean(cellfun(@(x) size(x, 1), temp)));
        wid = cat(1, wid, mean(cellfun(@(x) mean(x(:, 1)), temp)));
        temp = cat(1, temp{:});
        if isempty(temp)
            loc = cat(1, loc, {[]});
        else
            loc = cat(1, loc, {temp(:, 2)});
        end
    end
end

clust_stacks = cat(1, rsc.clust_stacks{2}, ee.clust_stacks{2});
r = nan(length(clust_stacks), 1);
for ii = 1:length(clust_stacks)
    if isempty(clust_stacks{ii})
        continue
    end
    temp = corr(clust_stacks{ii});
    idx = triu(ones(length(temp)), 1);
    r(ii) = mean(temp(~~idx));
end

% panel a - mean num place fields
g = iscue2(iscue2 ~= 0);
figure
subplot(1, 2, 1)
cdfplot(num_pf(iscue2 == 1));
hold on
cdfplot(num_pf(iscue2 == 2));

subplot(1, 2, 2)
boxplot(num_pf(iscue2 ~= 0), g);
ylim([1 2.5])

[~, p] = kstest2(num_pf(iscue2 == 1), num_pf(iscue2 == 2))
p = ranksum(num_pf(iscue2 == 1), num_pf(iscue2 == 2))

% panel b - mean place field width
g = iscue2(iscue2 ~= 0);
figure
subplot(1, 2, 1)
cdfplot(wid(iscue2 == 1));
hold on
cdfplot(wid(iscue2 == 2));

subplot(1, 2, 2)
boxplot(wid(iscue2 ~= 0), g);
ylim([5 35])

[~, p] = kstest2(wid(iscue2 == 1), wid(iscue2 == 2))
p = ranksum(wid(iscue2 == 1), wid(iscue2 == 2))

% panel c - fraction place cells
figure
boxplot(frac, iscue2)
ylim([0 1])
[p, ~, stats] = kruskalwallis(frac, iscue2)
multcompare(stats, 'ctype', 'bonferroni')

% panel d - place field locations
loc_cue = cat(1, loc{iscue2 == 1});
loc_traj = cat(1, loc{iscue2 == 2});

g = repelem((1:2)', [length(loc_cue); length(loc_traj)]);
loc_cat = [loc_cue; loc_traj];

g = gramm('x', loc_cat, 'color', g);
g.stat_bin('nbins', 20, 'geom', 'line', 'fill', 'transparent', 'normalization', 'probability');
g.axe_property('ylim', [0 .2]);
figure
g.draw();

% panel e - temporal correlation
g = iscue2(iscue2 ~= 0);
figure
subplot(1, 2, 1)
cdfplot(r(iscue2 == 1));
hold on
cdfplot(r(iscue2 == 2));

subplot(1, 2, 2)
boxplot(r(iscue2 ~= 0), g);
ylim([-.2 1])

[~, p] = kstest2(r(iscue2 == 1), r(iscue2 == 2))
p = ranksum(r(iscue2 == 1), r(iscue2 == 2))







