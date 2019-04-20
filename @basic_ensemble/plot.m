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
        deconv=fast_smooth(deconv(:,order),obj.ops.sig*obj.fs);
        ax(1)=subplot(5,1,1:3);
        imagesc('xdata',obj.ts,'cdata',deconv');
        colormap(get_colour('black'));
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
        for ii = 1:length(obj.clust_SCE)
            subset(ii)=plot(obj.ts,obj.clust_MUA(ii).MUA);
            hold on
            plot(obj.clust_SCE(ii).SCE.on, min(obj.clust_MUA(ii).MUA).*ones(length(obj.clust_SCE(ii).SCE.on),1),'^k');
            for i=1:length(obj.clust_SCE(ii).SCE.dur)
                plot(obj.clust_SCE(ii).SCE.on(i):.01:obj.clust_SCE(ii).SCE.on(i)+obj.clust_SCE(ii).SCE.dur(i), min(obj.clust_MUA(ii).MUA).*ones(length(obj.clust_SCE(ii).SCE.on(i):.01:obj.clust_SCE(ii).SCE.on(i)+obj.clust_SCE(ii).SCE.dur(i)),1), 'r-');
                text(obj.clust_SCE(ii).SCE.on(i)+obj.clust_SCE(ii).SCE.dur(i),min(obj.clust_MUA(ii).MUA),num2str(i));
            end
            ylabel('population mean dF/F');
            xlabel('time (sec)');
        end
        legend(subset,strsplit(num2str(1:length(obj.clust_SCE))));
        
        linkaxes(ax,'x');
        
    case 'lfp'
        if isempty(obj.lfp)
            error('need to load lfp first');
        end
        figure;
        deconv=fast_smooth(obj.deconv(:,order),obj.ops.sig*obj.fs);
        deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');
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
        
    case 'tree'
        if isempty(obj.tree)
            error('You need to build linkage tree first by running @basic_ensemble/hclust');
        end
        clist = 'brgycm';
        figure;
        ax_u=subplot(2,2,2);
        h = dendrogram(obj.tree,0,'reorder',order);
        set(h, 'color', 'k');
        for i = 1:length(obj.clust)
            idx = dendro_colours(obj.tree, obj.clust{i});
            set(h(idx), 'color', clist(mod(i,length(clist))+1));
        end
        axis square
        ax_l=subplot(2,2,3);
        h = dendrogram(obj.tree,0,'reorder',order(end:-1:1),'orientation','left');
        set(h, 'color', 'k');
        for i = 1:length(obj.clust)
            idx = dendro_colours(obj.tree, obj.clust{i});
            set(h(idx), 'color', clist(mod(i,length(clist))+1));
        end
        axis square
        ax_m=subplot(2,2,4);
        imagesc(obj.R(order,order));
        colormap jet
        axis square
        p = get(ax_m, 'Position');
        caxis([-.5 1])
        colorbar;
        set(ax_m, 'Position', p);
        
        linkaxes([ax_u ax_m], 'x');
        linkaxes([ax_l ax_m], 'y');
        
        
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
        
    case 'swr_window'
        if isempty(obj.swr_stack)
            error('need to run @basic_ensemble/swr_window first');
        end
        figure
        [~,idx] = max(obj.swr_stack);
        [~,idx] = sort(idx);
        imagesc('xdata',obj.swr_t,'cdata',zscore(fast_smooth(obj.swr_stack(:,idx),2))');
        xlim([min(obj.swr_t) max(obj.swr_t)]);
        ylim([1 size(obj.swr_stack,2)]);
        colormap jet
        c = colorbar;
        c.Label.String = 'normalized mean dF/F';
        title('Cortical Response Sorted');
        
        figure
        imagesc('xdata',obj.swr_t,'cdata',zscore(fast_smooth(obj.swr_stack(:,obj.order),2))');
        xlim([min(obj.swr_t) max(obj.swr_t)]);
        ylim([1 size(obj.swr_stack,2)]);
        colormap jet
        c = colorbar;
        c.Label.String = 'normalized mean dF/F';
        if strcmpi(obj.ops.order,'cluster')
            title('Cluster Sorted');
        elseif strcmpi(obj.ops.order,'pc')
            title('Place Field Sorted');
        else
        end
        xlabel('time from SWR peak (sec)');
        ylabel('sorted neuron no.');
        
    case 'silhouette'
        if strcmp(obj.ops.clust_method, 'shuffle')
            warning('Current clustering method set to ''shuffle'', changing method to ''silhouette''');
            obj.set_ops('clust_method','silhouette');
        end
        [~,s,ticks] = obj.silhouette_cluster;
        figure;
        plot(1:length(s), s, 'k');
        xticks(1:length(s));
        xticklabels(strsplit(num2str(ticks)));
        xlabel('number of clusters');
        ylabel('average silhouette width');
        
        
end