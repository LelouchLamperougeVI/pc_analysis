function h = single_cell_plot(analysis,k)
bins=50;
h=figure;
imagesc(analysis.psth{k});
set(gca,'xtick',0:bins/4:bins);
set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/4:0)));
colormap hot
ylabel('trials')
xlabel('distance (cm)')
colorbar