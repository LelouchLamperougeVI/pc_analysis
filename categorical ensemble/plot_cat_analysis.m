function plot_cat_analysis(analysis)

figure;
count=1;
for k=1:size(analysis.stack,3)
    if count>25
        count=1;
        figure;
    end
    subplot(5,5,count);
    imagesc(fast_smooth(analysis.stack(:,:,k)',1.5)');
    title(['n = ' num2str(k)]);
    colormap hot
    ylabel('trials')
    xlabel('time (arb.)')
    colorbar
    count=count+1;
end