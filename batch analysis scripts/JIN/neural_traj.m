fs = 1/ median(diff(tcs.tt));
speed = behavior.speed_raw;
epochs = abs(speed) > 0;
epochs = epochs - ~~movsum(epochs, round(fs * .25)) + epochs;

[U,S,V] = svd(zscore(fast_smooth(deconv, fs)));
traj = U(:, 1:2) * S(1:2, 1:2);

figure
% gscatter(traj(:, 1), traj(:, 2), epochs)
scatter(traj(:, 1), traj(:, 2), [], log(abs(speed)))