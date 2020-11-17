clear all
tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro');

test = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');
test.set_ops('e_size',5);
test.set_ops('clust_method','thres');
test.set_ops('sig', .2);
test.remove_mvt;
test.cluster;
deconv = test.twop.deconv;
deconv(any(isnan(deconv), 2), :) = [];
norm_deconv = (deconv - mean(deconv)) ./ std(deconv);

[~, l] = marchenkoPastur(norm_deconv);
[coeff, score, latent] = pca(norm_deconv);
sig_comps = latent > max(l);
pca_deconv = score(:, sig_comps) * coeff(:, sig_comps)';

[W, z, reconstructed, init] = rica_ensembles(pca_deconv, test.ensembles.clust);
z = (z - min(z)) ./ range(z);

% init = init(:, 1:7);
% md = rica(norm_deconv, size(init, 2), 'initialtransformweights', init, 'contrastfcn', 'logcosh');
% ev_thres = 0:.1:1;
% temp = [];
% for e = ev_thres
%     [W, z, reconstructed, init] = rica_ensembles(norm_deconv, test.ensembles.clust, e);
%     temp = cat(2, temp, diag(W' * init));
% end

figure
ax(1) = subplot(3,1,1);
smthnd = fast_smooth(norm_deconv(:, test.analysis.order), 5)';
imagesc(smthnd)
% caxis([prctile(smthnd(:), 0.01) prctile(smthnd(:), 99.9)])
ax(2) = subplot(3,1,2);
smthnd = fast_smooth(reconstructed(:, test.analysis.order), 5)';
imagesc(smthnd)
% caxis([prctile(smthnd(:), 0.01) prctile(smthnd(:), 99.9)])
ax(3) = subplot(3,1,3);
plot(z + (1:size(z,2)));
linkaxes(ax, 'x')

figure
subplot(2,1,1)
plot(init(test.analysis.order, 1:length(test.ensembles.clust)));
legend
subplot(2,1,2)
plot(W(test.analysis.order, :));
legend

diag(corr(init(test.analysis.order, 1:length(test.ensembles.clust)), W))
diag(W' * init)

%%%
beh_dec = test.analysis.original_deconv;
beh_dec = (beh_dec - min(beh_dec)) ./ range(beh_dec);
score = beh_dec * W;
trial_bins = discretize(1:size(score,1), test.analysis.behavior.trials_ts);
trial_bins(isnan(trial_bins)) = 0;
pos_bins = discretize(test.analysis.behavior.unit_pos, length(test.analysis.Pi));

psth = arrayfun(@(x) accumarray([pos_bins' trial_bins' + 1], score(:,x), [length(unique(pos_bins)) length(unique(trial_bins))], @mean), 1:size(score,2), 'uniformoutput', false);
psth = cellfun(@(x) x(:, 2:end), psth, 'uniformoutput', false);

for ii = 1:length(psth)
    if mod(ii, 9) == 1
        figure
    end
    subplot(3, 3, mod(ii - 1, 9) + 1);
    imagesc(psth{ii}');
    colorbar
    colormap hot
end


%%%
% test = rand(1000, 3);
% md = rica(test, 2);
% z = transform(md, test);


