clear all

animal = 'EE006';
date = '2019_06_05';

lfp = ensemble(fullfile('/mnt/storage/rrr_magnum/M2/', animal, date, [date '_3.abf']));
lfp.set_ops('e_size',5);
lfp.set_ops('clust_method','thres');
lfp.set_ops('sig', .2);
lfp.remove_mvt;
lfp.cluster;
lfp.set_ops('order','cluster')
lfp.detect_sce;
lfp.hPICA;

targets = [4, 6];
idx = [];
for ii = 1:length(targets)
    idx = cat(1, idx, lfp.ensembles.clust{targets(ii)}');
end
idx = cat(1, idx, randsample(setxor(1:size(lfp.twop.deconv, 2), cell2mat(lfp.ensembles.clust)), 20)');

deconv = lfp.twop.deconv(:, idx);
deconv = fast_smooth(deconv, lfp.ops.sig * lfp.twop.fs);
deconv = (deconv - nanmin(deconv)) ./ range(deconv);

R = corr(deconv(~any(isnan(deconv), 2), :));

D = 1 - R;
D = pdist(D);
tree = linkage(D, 'average');
order = optimalleaforder(tree, D);

clear ax
figure
ax(2, 1) = subplot(2, 2, 3);
imagesc(-deconv(:, order)')
colormap(ax(2, 1), bone)
xlim([1.1e4, 1.6e4])
ax(1, 2) = subplot(2, 2, 2);
dendrogram(tree, 0, 'reorder', order);
axis square
ax(2, 2) = subplot(2, 2, 4);
imagesc(R(order, order));
colormap(ax(2, 2), jet)
axis square
linkaxes(ax(2, :), 'y')
linkaxes(ax(:, 2), 'x')

clear ax
figure
ax(1) = subplot(2, 1, 1);
imagesc(-deconv(:, order)')
colormap bone
ax(2) = subplot(2, 1, 2);
plot(lfp.hiepi.z(:, targets))
linkaxes(ax, 'x')
xlim([1.1e4, 1.6e4])


figure
k = length(targets);
img = [];
for ii = 1:k
    subplot(2, k, ii);
    clust = lfp.analysis.stack(:, lfp.ensembles.clust{targets(ii)});
    [~, idx] = max(clust);
    [~, idx] = sort(idx);
    imagesc(-clust(:, idx)')
    img{ii} = -clust(:, idx)';
    colormap bone
    
    subplot(2, k, ii + k);
    clust = lfp.hiepi.psth{targets(ii)};
    clust = mean(clust, 2);
    clust = fast_smooth(clust, 5 / 150 * 50);
    plot(clust);
end


