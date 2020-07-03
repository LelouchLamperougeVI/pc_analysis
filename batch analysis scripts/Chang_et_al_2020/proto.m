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