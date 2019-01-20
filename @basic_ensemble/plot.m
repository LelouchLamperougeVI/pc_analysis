function plot(obj,type)
%Plot results from detection

figure;
if ~isempty(obj.order)
    order=obj.order;
    lab='sorted neuron no.';
else
    order=1:size(obj.deconv,2);
    lab='neuron no.';
end

switch lower(type)
    case 'sce'
        deconv=fast_smooth(obj.deconv(:,order),obj.ops.sig);
        ax(1)=subplot(4,1,1:3);
        imagesc('xdata',obj.ts,'cdata',deconv');
        ylim([1 size(deconv,2)]);
        ylabel(lab);

        ax(2)=subplot(4,1,4);
        plot(obj.ts,obj.MUA);
        hold on
        plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
        for i=1:length(obj.SCE.dur)
            plot(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i), min(obj.MUA).*ones(length(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i)),1), 'r-');
        end
        ylabel('population mean dF/F');
        xlabel('time (sec)');

        linkaxes(ax,'x');
        
    case 'corr'
        if isempty(obj.R)
            error('no correlation matrix available');
        end
        subplot(1,2,1)
        imagesc(obj.R(order,order));
        colormap jet
        axis square
        caxis([-.2 1]);
        c=colorbar;
        c.Label.String='Corr. Coef.';
        xlabel(lab);
        ylabel(lab);
        title('SCE');
        
        subplot(1,2,2)
        imagesc(obj.null_R(order,order));
        colormap jet
        axis square
        caxis([-.2 1]);
        c=colorbar;
        c.Label.String='Corr. Coef.';
        xlabel(lab);
        ylabel(lab);
        title('outside detection');
        
end