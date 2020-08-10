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
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest1.ops.sig * rest1.twop.fs);
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
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest1.ops.sig * run.twop.fs);
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
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest2.ops.sig * rest2.twop.fs);
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
rest2.set_ops('e_size',5);
rest2.set_ops('clust_method','thres');
rest2.cluster;
rest2.topography;

rest2.plot('clust_topo')