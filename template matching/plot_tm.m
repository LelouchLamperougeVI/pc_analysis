function plot_tm(analysis,deconv)
black=[linspace(1,0,64)' linspace(1,0,64)' linspace(1,0,64)'];

order=analysis.stack;
[~,order]=max(order);
[~,order]=sort(order);

figure;
imagesc(analysis.template(:,order)');
title('template');
colormap(black);

figure;
ax1=subplot(4,1,1:2);
imagesc(fast_smooth(zscore(deconv(:,order)),2)');
colormap(ax1,black);
c=colorbar; c.Label.String='dF/F';

ax2=subplot(4,1,3);
imagesc(analysis.C');
colormap(ax2,hot);
c=colorbar; c.Label.String='Spearman''s Rho';

ax3=subplot(4,1,4);
imagesc(analysis.C_pval');
colormap(ax3,bone);
c=colorbar; c.Label.String='prctile';
[x,y]=find(analysis.C_pval>.95);
hold on
plot(x,y,'r*');

linkaxes([ax1 ax2 ax3], 'x');