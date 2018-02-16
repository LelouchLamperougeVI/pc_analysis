function plot_analysis(analysis,plotFlag)

bins=50;

if nargin==1
    plotFlag=[1 1];
end

if plotFlag(1)
    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
    ordered=ordered(any(analysis.pc_list'==ordered));
    
    figure;
    count=1;
    for k=ordered
        if 1
            if count>25
                count=1;
                figure;
            end
            subplot(5,5,count);
            imagesc(log(analysis.psth{k}));
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
end

if plotFlag(2)
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
end