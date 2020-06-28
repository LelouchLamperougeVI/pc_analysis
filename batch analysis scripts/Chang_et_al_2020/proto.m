clear all

root = '/mnt/storage/rrr_magnum/M2';
animal = 'EE001';
date = '2018_12_28';
session = '3';

% chan = [1 2 3 5 NaN 6 NaN];
chan = [1 2 3 6 NaN 5 NaN];

ass = ensemble(fullfile(root, animal, date, session, [date '_' session '.abf']), 'channels', chan);

%%
ass.set_ops('e_size',10);
ass.set_ops('clust_method','thres');
ass.cluster;
ass.set_ops('order','pc');

%% plot rest-run
figure
subplot(1, 2, 1);
deconv = ass.analysis.deconv(:, ass.ensembles.clust_order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(deconv');
colormap hot
subplot(1, 2, 2);
deconv = ass.twop.deconv(:, ass.ensembles.clust_order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(deconv');
colormap hot


%% make batch analysis
clear all

root = '/mnt/storage/rrr_magnum/M2';
animal = 'EE001';
dates = {'2018_12_18', '2018_12_20', '2018_12_28', '2019_01_07'}';

for ii = 1:length(dates)
    lfp = lfp(fullfile(root, animal, dates{ii}, [dates{ii} '_2.abf']));
    lfp.perform_analysis;
    lfp.save('analysis');
    lfp.save('chan');
end