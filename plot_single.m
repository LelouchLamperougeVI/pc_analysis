function h=plot_single(analysis,k)

bins=length(analysis.Pi);

for i=1:length(k)
    h=figure;
    imagesc((analysis.psth{k(i)}));
    set(gca,'xtick',0:bins/4:bins);
    set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/4:0)));
    title(['n = ' num2str(k(i))]);
    colormap hot
    ylabel('trials')
    xlabel('distance (cm)')
    colorbar
end