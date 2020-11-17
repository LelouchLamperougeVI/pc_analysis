clear all

test = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');
test.set_ops('e_size',5);
test.set_ops('clust_method','thres');
test.set_ops('sig', .2);
test.remove_mvt;
test.cluster;
deconv = test.twop.deconv;
deconv(any(isnan(deconv), 2), :) = [];

%% choose test neurons
norm_deconv = (deconv - mean(deconv)) ./ std(deconv);
% figure
% imagesc(fast_smooth(norm_deconv(:, test.analysis.order), 5)');
% test_neurs = test.analysis.order([109 110 111 29 30 31 95 96 97]);
test_neurs = test.analysis.order;
norm_deconv = norm_deconv(:, test_neurs);

%% perform pca
[coeff, score, latent] = pca(norm_deconv);

[~, l] = marchenkoPastur(norm_deconv);
sig_comp = latent > max(l);
pca_deconv = score(:, sig_comp) * coeff(:, sig_comp)';

%% perform rica
init = zeros(size(pca_deconv, 2), length(test.ensembles.clust)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(test.ensembles.clust)
    init(test.ensembles.clust{ii}, ii) = 1;
end
init = init ./ vecnorm(init);
init = init(test.analysis.order, :);

pca_deconv = (pca_deconv - mean(pca_deconv)) ./ std(pca_deconv);
md = rica(pca_deconv, 7, 'InitialTransformWeights', init);
rica_deconv = pca_deconv * (md.TransformWeights * md.TransformWeights');

diag(md.TransformWeights' * init)

figure
ax(1) = subplot(3,1,1);
imagesc(norm_deconv')
ax(2) = subplot(3,1,2);
imagesc(pca_deconv')
ax(3) = subplot(3,1,3);
imagesc(rica_deconv')
linkaxes(ax, 'x')

figure
plot(md.TransformWeights)

figure
ax(1) = subplot(3,1,1);
imagesc(norm_deconv')
ax(2) = subplot(3,1,2);
temp = norm_deconv * coeff(:, sig_comp);
plot((temp - min(temp)) ./ range(temp) + (1:sum(sig_comp)));
ax(3) = subplot(3,1,3);
temp = pca_deconv * md.TransformWeights;
plot((temp - min(temp)) ./ range(temp) + (1:7));
linkaxes(ax, 'x')

figure
subplot(1,3,1)
imagesc(corr(norm_deconv))
axis square
subplot(1,3,2)
imagesc(corr(pca_deconv))
axis square
subplot(1,3,3)
imagesc(corr(rica_deconv))
axis square