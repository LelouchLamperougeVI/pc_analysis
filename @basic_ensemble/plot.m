function plot(obj,type,varargin)
%Plot results from detection

if ~isempty(obj.order)
    order=obj.order;
    lab='sorted neuron no.';
else
    order=1:size(obj.deconv,2);
    lab='neuron no.';
end

switch lower(type)
    case 'sce'
        figure;
        deconv=(obj.deconv - mean(obj.deconv,'omitnan'))./std(obj.deconv,'omitnan');
        deconv=fast_smooth(deconv(:,order),obj.ops.sig);
        ax(1)=subplot(4,1,1:3);
        imagesc('xdata',obj.ts,'cdata',deconv');
        colormap(get_colour('black'));
        ylim([1 size(deconv,2)]);
        ylabel(lab);
        
        ax(2)=subplot(4,1,4);
        plot(obj.ts,obj.MUA);
        hold on
        plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
        for i=1:length(obj.SCE.dur)
            plot(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i), min(obj.MUA).*ones(length(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i)),1), 'r-');
            text(obj.SCE.on(i)+obj.SCE.dur(i),min(obj.MUA),num2str(i));
        end
        ylabel('population mean dF/F');
        xlabel('time (sec)');
        
        linkaxes(ax,'x');
        
    case 'lfp'
        if isempty(obj.lfp)
            error('need to load lfp first');
        end
        figure;
        deconv=fast_smooth(obj.deconv(:,order),obj.ops.sig);
%         deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');
        ax(1)=subplot(5,1,1:3);
        imagesc('xdata',obj.ts,'cdata',deconv');
%         colormap(get_colour('black'));
        ylim([1 size(deconv,2)]);
        ylabel(lab);
        
        ax(2)=subplot(5,1,4);
        plot(obj.ts,obj.MUA);
        hold on
        plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
        for i=1:length(obj.SCE.dur)
            plot(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i), min(obj.MUA).*ones(length(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i)),1), 'r-');
            text(obj.SCE.on(i)+obj.SCE.dur(i),min(obj.MUA),num2str(i));
        end
        ylabel('population mean dF/F');
        xlabel('time (sec)');
        
        ax(3)=subplot(5,1,5);
        plot(obj.lfp.t, obj.lfp.lfp);
        xlabel('time (sec)');
        ylabel('lfp');
        
        linkaxes(ax,'x');
        
    case 'corr'
        if isempty(obj.R)
            error('no correlation matrix available');
        end
        figure;
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
        
    case 'clust_corr'
        if isempty(obj.clust)
            error('data not clustered');
        end
        k=3;
        count=1;
        while(count<=length(obj.clust))
            if ~mod(count-1,k^2)
                figure;
            end
            subplot(k,k,mod(count-1,k^2)+1);
            R=obj.R(order,order);
            [~,idx]=intersect(order, setxor(order,obj.clust{count}));
            R(idx,idx)=0;
            imagesc(R);
            colormap jet
            axis square
            caxis([0 .6]);
            c=colorbar;
            c.Label.String='Corr. Coef.';
            xlabel(lab);
            ylabel(lab);
            title(['clust ' num2str(count)]);
            count=count+1;
        end
        
    case 'pc'
        if isempty(obj.clust)
            error('data not clustered');
        end
        if nargin<3
            error('please give cluster number');
        end
        plot_analysis(obj.lfp.analysis,[1 0 0],obj.clust{varargin{1}});
        
    case 'spec'
        if isempty(obj.spec)
            error('need to compute spectrum first')
        end
        figure;
        subplot(1,2,1);
        spec=(obj.spec.spectrum_on-mean(obj.spec.norm,2))./std(obj.spec.norm,0,2);
        spec=imgaussfilt(spec,1);
        imagesc('xdata',obj.spec.t,'ydata',obj.spec.f,'cdata',spec);
        colormap jet
        caxis([0 max(spec(:))]);
        xlim([obj.spec.t(1) obj.spec.t(end)]);
        ylim([obj.spec.f(1) obj.spec.f(end)]);
        c=colorbar;
        c.Label.String='zscore power';
        xlabel('time from SCE onset (sec)');
        ylabel('frequency');
        title('SCE onset')
        
        subplot(1,2,2);
        spec=(obj.spec.spectrum_peak-mean(obj.spec.norm,2))./std(obj.spec.norm,0,2);
        spec=imgaussfilt(spec,1);
        imagesc('xdata',obj.spec.t,'ydata',obj.spec.f,'cdata',spec);
        colormap jet
        caxis([0 max(spec(:))]);
        xlim([obj.spec.t(1) obj.spec.t(end)]);
        ylim([obj.spec.f(1) obj.spec.f(end)]);
        c=colorbar;
        c.Label.String='zscore power';
        xlabel('time from SCE peak (sec)');
        ylabel('frequency');
        title('SCE peak');
        
end