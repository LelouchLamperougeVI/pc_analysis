function plot_tm(analysis,deconv,cf)
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

if nargin>2
    C=analysis.C_pval(:,cf)>.95;
    idx=size(analysis.template,1)/cf;
    idx=conv(C,ones(round(idx*2)+1,1),'same')>0;
    
    figure;
    deconv(~idx,:)=max(deconv(:));
    imagesc(deconv(:,order)');
    colormap(black);
    title(['Significant epochs at CF = ' num2str(cf)]);
end