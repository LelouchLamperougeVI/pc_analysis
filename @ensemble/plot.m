function plot(obj,type,varargin)
%Plot results from detection

if ~isempty(obj.order)
    order=obj.order;
    lab='sorted neuron no.';
else
    order=1:size(obj.twop.deconv,2);
    lab='neuron no.';
end

switch lower(type)
    case 'sce'
        wbins = 100;
        figure;
        %         deconv=(obj.twop.deconv - mean(obj.twop.deconv,'omitnan'))./std(obj.twop.deconv,'omitnan');
        %         deconv=fast_smooth(deconv(:,order),obj.ops.sig*obj.twop.fs);
        deconv=fast_smooth(obj.twop.deconv(:,order),obj.ops.sig*obj.twop.fs);
        %         deconv=(deconv - mean(deconv,'omitnan'))./std(deconv,'omitnan');
        deconv = (deconv - min(deconv,[],'omitnan')) ./ range(deconv);
        idx = 1:3*wbins;
        idx( ~mod(1:3*wbins, wbins) ) = [];
        ax(1)=subplot(5,wbins,idx);
        imagesc('xdata',obj.twop.ts,'cdata',deconv');
        %         colormap(get_colour('black'));
        colormap hot
        ylim([1 size(deconv,2)]);
        ylabel(lab);
        
        idx = 3*wbins:4*wbins;
        idx( ~mod(3*wbins:4*wbins, wbins) ) = [];
        ax(2)=subplot(5,wbins,idx);
        plot(obj.twop.ts,obj.MUA);
        hold on
        plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
        for i=1:length(obj.SCE.dur)
            plot(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i), min(obj.MUA).*ones(length(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i)),1), 'r-');
            text(obj.SCE.on(i)+obj.SCE.dur(i),min(obj.MUA),num2str(i));
        end
        ylabel('population mean dF/F');
        xlabel('time (sec)');
        
        idx = 4*wbins:5*wbins;
        idx( ~mod(4*wbins:5*wbins, wbins) ) = [];
        ax(3)=subplot(5,wbins,idx);
        subset = gobjects(length(obj.clust_SCE), 1);
        for ii = 1:length(obj.clust_SCE)
            subset(ii)=plot(obj.twop.ts,obj.clust_MUA(ii).MUA, 'color', obj.colours(ii,:));
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
        
        idx = 1:3*wbins;
        idx( ~~mod(1:3*wbins, wbins) ) = [];
        ax(4) = subplot(5,wbins, idx);
        if strcmpi(obj.ops.order,'cluster')
            hold on
            for ii = 1:length(obj.clust)
                plot(zeros(1, length(obj.clust{ii})), arrayfun(@(x) find(obj.order == x), obj.clust{ii}), 'linewidth',8, 'color',obj.colours(ii,:));
            end
        end
        linkaxes([ax(1) ax(4)],'y');
        set(ax(4), 'visible','off');
        
    case 'lfp'
        if isempty(obj.lfp.lfp)
            error('need to load lfp first');
        end
        figure;
        deconv=fast_smooth(obj.twop.deconv(:,order),obj.ops.sig*obj.twop.fs);
        deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');
        ax(1)=subplot(5,1,1:3);
        imagesc('xdata',obj.twop.ts,'cdata',deconv');
        %         colormap(get_colour('black'));
        ylim([1 size(deconv,2)]);
        ylabel(lab);
        
        ax(2)=subplot(5,1,4);
        plot(obj.twop.ts,obj.MUA);
        hold on
        plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
        for i=1:length(obj.SCE.dur)
            plot(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i), min(obj.MUA).*ones(length(obj.SCE.on(i):.01:obj.SCE.on(i)+obj.SCE.dur(i)),1), 'r-');
            text(obj.SCE.on(i)+obj.SCE.dur(i),min(obj.MUA),num2str(i));
        end
        ylabel('population mean dF/F');
        xlabel('time (sec)');
        
        ax(3)=subplot(5,1,5);
        plot(obj.lfp.ts, obj.lfp.lfp);
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
            title(['clust ' num2str(count)], 'color',obj.colours(count,:));
            count=count+1;
        end
        
    case 'tree'
        if isempty(obj.tree)
            error('You need to build linkage tree first by running @basic_ensemble/hclust');
        end
        figure;
        ax_u=subplot(2,2,2);
        h = dendrogram(obj.tree,0,'reorder',order);
        set(h, 'color', 'k');
        for i = 1:length(obj.clust)
            idx = dendro_colours(obj.tree, obj.clust{i});
            set(h(idx), 'color', obj.colours(i,:));
        end
        axis square
        ax_l=subplot(2,2,3);
        h = dendrogram(obj.tree,0,'reorder',order(end:-1:1),'orientation','left');
        set(h, 'color', 'k');
        for i = 1:length(obj.clust)
            idx = dendro_colours(obj.tree, obj.clust{i});
            set(h(idx), 'color', obj.colours(i,:));
        end
        axis square
        ax_m=subplot(2,2,4);
        imagesc(obj.R(order,order));
        
        rbmap(ax_m, 'colorbar',true, 'caxis', [-.5 1], 'normalize',true, 'interp', 129);
        axis square
        
        linkaxes([ax_u ax_m], 'x');
        linkaxes([ax_l ax_m], 'y');
        
    case 'pc'
        if isempty(obj.clust)
            error('data not clustered');
        end
        if nargin<3
            error('please give cluster number');
        end
        plot_analysis(obj.analysis,[1 1 0],obj.clust{varargin{1}});
        
    case 'spec'
        figure;
        subplot(1,2,1);
        spec=(obj.lfp.spec.spectrum_on-mean(obj.lfp.spec.norm,2))./std(obj.lfp.spec.norm,0,2);
        spec=imgaussfilt(spec,1);
        imagesc('xdata',obj.lfp.spec.t,'ydata',obj.lfp.spec.f,'cdata',spec);
        colormap jet
        caxis([0 max(spec(:))]);
        xlim([obj.lfp.spec.t(1) obj.lfp.spec.t(end)]);
        ylim([obj.lfp.spec.f(1) obj.lfp.spec.f(end)]);
        c=colorbar;
        c.Label.String='zscore power';
        xlabel('time from SCE onset (sec)');
        ylabel('frequency');
        title('SCE onset')
        
        subplot(1,2,2);
        spec=(obj.lfp.spec.spectrum_peak-mean(obj.lfp.spec.norm,2))./std(obj.lfp.spec.norm,0,2);
        spec=imgaussfilt(spec,1);
        imagesc('xdata',obj.lfp.spec.t,'ydata',obj.lfp.spec.f,'cdata',spec);
        colormap jet
        caxis([0 max(spec(:))]);
        xlim([obj.lfp.spec.t(1) obj.lfp.spec.t(end)]);
        ylim([obj.lfp.spec.f(1) obj.lfp.spec.f(end)]);
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
            hold on
            for ii = 1:length(obj.clust)
                plot(zeros(1, length(obj.clust{ii})), arrayfun(@(x) find(obj.order == x), obj.clust{ii}), 'linewidth',5, 'color',obj.colours(ii,:));
            end
        elseif strcmpi(obj.ops.order,'pc')
            title('Place Field Sorted');
        else
        end
        xlabel('time from SWR peak (sec)');
        ylabel('sorted neuron no.');
        
        k=5;
        h = gobjects(length(obj.clust), 1);
        for ii=1:length(obj.clust)
            if ~mod(ii-1, k^2)
                figure;
            end
            h(ii) = subplot(k,k, mod(ii-1, k^2)+1);
            temp = mean(obj.swr_all(:, obj.clust{ii}, :), 2);
            err = sem(temp, 3);
            mu = mean(temp, 3);
            errorshade(obj.swr_t, mu, err, 'colour', obj.colours(ii,:), 'target',h(ii));
            temp = mean(obj.swr_null_all(:, obj.clust{ii}, :), 2);
            err = sem(temp, 3);
            mu = mean(temp, 3);
            errorshade(obj.swr_t, mu, err, 'colour', 'k', 'target',h(ii));
            title(['Clust ' num2str(ii)]);
        end
        linkaxes(h, 'y');
        %         labels = arrayfun(@(x) ['clust ' num2str(x)], 1:length(obj.clust), 'uniformoutput',false);
        %         legend(labels);
        
        for ii=1:length(obj.clust)
            if ~mod(ii-1, k^2)
                figure;
            end
            h(ii) = subplot(k,k, mod(ii-1, k^2)+1);
            temp = mean(obj.swr_all(:, obj.clust{ii}, :), 2);
            %             imagesc(squeeze(temp)');
            imagesc(fast_smooth(squeeze(temp),1)');
            title(['Clust ' num2str(ii)]);
        end
        
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
        
    case 'clust_topo'
        if isempty(varargin)
            clusts = 1:length(obj.clust);
        else
            clusts = varargin{1};
        end
        k=2;
        for plane = 1:size(obj.topo.mimg, 3)
            if ~mod(plane-1, k^2)
                figure;
            end
            subplot(k, k, mod(plane-1, k^2) + 1);
            
            mimg = obj.topo.mimg(:, :, plane);
            mimg = (mimg - min(mimg(:))) ./ range(mimg(:));
            imshow(mimg)
            hold on
            
            mask = zeros(size(obj.topo.maskNeurons, 1), size(obj.topo.maskNeurons, 2), 3);
            for ii = clusts
                idx = ismember(obj.topo.maskNeurons(:, :, plane), obj.clust{ii});
                idx = repmat(idx, 1, 1, 3);
                mask = mask + double(idx) .* permute(obj.colours(ii, :), [1 3 2]);
            end
            
            h = imshow(mask);
            set(h, 'alphadata', .5);
            title(obj.twop.planes.plane_names{obj.twop.planes.planes(plane)});
        end
        
    case 'clust_topo_stack'
        if isempty(varargin)
            clusts = 1:length(obj.clust);
        else
            clusts = varargin{1};
        end
        
        figure
        hold on
        x = linspace(0, obj.topo.FOV(1), size(obj.topo.mimg, 1));
        y = linspace(0, obj.topo.FOV(2), size(obj.topo.mimg, 2));
        sections = -unique(obj.twop.planes.depth);
        [x, y, z] = meshgrid(y, x, sections);
        slice(x, y, z, obj.topo.mimg, [], [], sections)
        axis image
        colormap gray
        shading interp
        alpha(.3)
        
        for ii = clusts
            plot3(obj.topo.centroid(1, obj.clust{ii}), obj.topo.centroid(2, obj.clust{ii}), -obj.twop.planes.depth(obj.clust{ii}), '.', 'color', obj.colours(ii, :), 'markersize', 30);
        end
        
        xlabel('galvo')
        ylabel('resonant')
        zlabel('piezo')
        
    case 'clust_topo_stats'
        figure
        subplot(2, 2, 1);
        cdfplot(cell2mat(obj.topo.clust.zsilhouette));
        xlim([-1 1]);
        xlabel('silhouette score')
        ylabel('cum. freq.')
        xline(0);
        
        subplot(2, 2, 2);
        violin(obj.topo.clust.zsilhouette, 'scatter');
        xlabel('ensemble');
        ylabel('silhouette score')
        yline(0);
        
        subplot(2, 2, 3);
        d = cellfun(@(x) squareform(obj.topo.zdistances(x, x)), obj.clust, 'uniformoutput', false);
        violin(d, 'box');
        xlabel('ensemble');
        ylabel('distance (\mum)')
        d = squareform(obj.topo.zdistances);
        yline(median(d));
        title('line: sample median')
        
        subplot(2, 2, 4);
        depth = obj.twop.planes.depth(cell2mat(obj.clust));
        s = cell2mat(obj.topo.clust.zsilhouette);
        s = arrayfun(@(x) s(depth == x), unique(depth), 'uniformoutput', false);
        violin(s, 'labels', strsplit(num2str(unique(depth))), 'box');
        ylim([-1 1]);
        xlabel('plane depth');
        ylabel('silhouette score')
        yline(0);
        
    otherwise
        plot@LFP(obj, type);
        
end