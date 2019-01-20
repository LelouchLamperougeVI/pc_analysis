function plot_analysis(analysis,plotFlag,pc_list)
% [PSTHs    stack/cMatrix   SI/pc_width]
if nargin<2
    plotFlag=[1 1 1];
end

bins=length(analysis.Pi);
if nargin < 3
    pc_list=analysis.pc_list;
end

if plotFlag(1)
    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
    ordered=ordered(any(pc_list'==ordered));
    
    figure;
    count=1;
    for k=ordered
        if 1
            if count>25
                count=1;
                figure;
            end
            subplot(5,5,count);
%             imagesc(log(analysis.psth{k}));
            imagesc((analysis.psth{k}));
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
    stack=analysis.stack(:,pc_list);

%     sd=4.5/analysis.vr_length*bins;
%     stack=analysis.raw_psth(:,:,pc_list);
%     stack=squeeze(mean(stack));
%     stack=fast_smooth(stack,sd);
%     stack=(stack-min(stack))./range(stack);
    
%     stack=arrayfun(@(x) mean(fast_smooth(analysis.raw_psth(:,:,x),sd,2)),pc_list,'uniformoutput',false);
%     stack=cell2mat(stack);
%     stack=reshape(stack,bins,length(pc_list));
%     stack=(stack-min(stack))./range(stack);
    
    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
%     ordered=ordered(any(pc_list'==ordered));
    stack=stack(:,ordered)';

    figure;
    imagesc(stack);
    set(gca,'xtick',0:bins/5:bins);
    set(gca,'xticklabel',strsplit(num2str(-analysis.vr_length:analysis.vr_length/5:0)));
    xlabel('position (cm)');
    ylabel('ordered neuron no.');
%     colormap jet; 
    red=[ones(64,1) linspace(1,0,64)' linspace(1,0,64)'];
    blue=red(:,[2 3 1]);
    black=[linspace(1,0,64)' linspace(1,0,64)' linspace(1,0,64)'];
    colormap(black); %red
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
%     colormap jet
    colormap(red);
    axis square
end

if plotFlag(3)
    pc_width=vertcat(analysis.width{:});
    pc_width=pc_width.*analysis.vr_length./length(analysis.Pi);
    
    figure;
    subplot(1,3,1);
    [f,x]=ecdf(analysis.SI);
    plot(x,f);
    xlabel('proficiency')
    ylabel('cumm. prob.')
    title('Uncertainty Coefficient');
    subplot(1,3,2);
    [f,x]=ecdf(analysis.sparsity);
    plot(x,f);
    xlabel('sparsity')
    ylabel('cumm. prob.')
    title('Sparsity');
    subplot(1,3,3);
    [f,x]=ecdf(pc_width(:,1));
    plot(x,f);
    xlabel('place fields width (cm)')
    ylabel('cumm. prob.')
    title('Place Fields Width');
end
