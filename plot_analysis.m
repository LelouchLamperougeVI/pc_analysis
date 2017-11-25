function plot_analysis(analysis)

bins=50;
figure;
count=1;
for k=analysis.pc_list
    if 1
        if count>25
            count=1;
            figure;
        end
        subplot(5,5,count);
        imagesc(analysis.psth{k});
        set(gca,'xtick',0:bins/4:bins);
        set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/4:0)));
        title(['n = ' num2str(k)]);
        colormap hot
        ylabel('trials')
        xlabel('distance (cm)')
        colorbar
    end
    count=count+1;
end

stack=analysis.raw_stack(:,analysis.pc_list);

stack=(stack-repmat(min(stack),bins,1));
stack=stack./repmat(max(stack),bins,1);

[~,idx]=max(stack);
[~,ordered]=sort(idx);
stack=stack(:,ordered)';
figure;
imagesc(stack);
set(gca,'xtick',0:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/5:0)));
xlabel('position (cm)');
ylabel('ordered neuron no.');
colormap jet;
c=colorbar; c.Label.String='Norm. Mean dF/F';


qMatrix=corr(stack);

figure;
imagesc(qMatrix);
set(gca,'xtick',0:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/5:0)));
xlabel('position (cm)');
set(gca,'ytick',0:bins/5:bins);
set(gca,'yticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/5:0)));
ylabel('position (cm)');
c=colorbar; c.Label.String='corr. coef.';
colormap jet
axis square