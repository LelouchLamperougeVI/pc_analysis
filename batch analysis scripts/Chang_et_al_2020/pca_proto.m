clear all

test = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');
test.set_ops('e_size',5);
test.set_ops('clust_method','thres');
test.set_ops('sig', .2);
test.remove_mvt;
test.cluster;
order = test.ensembles.clust_order;
rmd = rica_ensembles(test.twop.deconv, test.ensembles.clust, 'denoise', true);
z = (rmd.z - min(rmd.z)) ./ range(rmd.z);

figure
ax(1) = subplot(3,1,1);
smthnd = fast_smooth(rmd.X(:, order), 5)';
imagesc(smthnd)
ax(2) = subplot(3,1,2);
smthnd = fast_smooth(rmd.reconst(:, order), 5)';
imagesc(smthnd)
ax(3) = subplot(3,1,3);
plot(z + (1:size(z, 2)));
linkaxes(ax, 'x')

figure
subplot(2,1,1)
plot(rmd.init(order, :));
legend
subplot(2,1,2)
plot(rmd.W(order, :));
legend

%%%
beh_dec = test.analysis.original_deconv;
beh_dec = (beh_dec - min(beh_dec)) ./ range(beh_dec);
score = beh_dec * rmd.W;
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
