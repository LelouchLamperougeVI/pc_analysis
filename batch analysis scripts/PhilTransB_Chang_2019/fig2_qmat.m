%% Run
figure
ax(1) = subplot(7,1,1:5);
imagesc('xdata',lfp.ts_2p,'cdata', fast_smooth(zscore(lfp.deconv(:,get_order(analysis))),30)');
colormap hot
ylim([1 size(lfp.deconv,2)]);
p = get(ax(1), 'Position');
caxis([0 2]);
colorbar;
set(ax(1), 'Position', p);

ax(2) = subplot(7,1,6);
plot(lfp.ts_2p, behavior.unit_pos);

ax(3) = subplot(7,1,7);
plot(lfp.ts_2p, behavior.unit_vel);

linkaxes(ax, 'x');

xlim([120 180])

%% Rest
figure
ax(1) = subplot(7,1,1:5);
imagesc('xdata',lfp.ts_2p,'cdata', fast_smooth(zscore(lfp.deconv(:,get_order(analysis))),30)');
% imagesc('xdata',lfp.ts_2p,'cdata', fast_smooth(zscore(lfp.deconv(:,ass.order)),30)');
colormap hot
ylim([1 size(lfp.deconv,2)]);
p = get(ax(1), 'Position');
caxis([0 2]);
colorbar;
set(ax(1), 'Position', p);

ax(2) = subplot(7,1,6);
% plot(lfp.ts_2p, behavior.unit_pos);

ax(3) = subplot(7,1,7);
plot(lfp.ts_2p, lfp.behavior.speed_raw_noSmooth);

linkaxes(ax, 'x');

% xlim([180 240]) %rest1
xlim([550 610]) %rest2
