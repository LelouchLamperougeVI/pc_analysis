function plot(obj,option)
% Options:
%   'spectrum'
%   'bands'
%   'deconv_spec'
%   'channels"

switch option
    case 'spectrum'
        if isempty(obj.lfp.spec.spectrum)
            warning('not spectrum was generated; attempting to generate spectrum');
            obj.spectrum;
        end
        figure;
        plot_spec;
        
    case 'bands'
        k=5;
        if isempty(obj.lfp.theta)
            error('must filter bands first');
        end
        %         if ~isempty(obj.camera.cam.mvt)
        %             k=k+1;
        %         end
        figure;
        ax1=subplot(k,1,1);
        plot(obj.lfp.ts,obj.lfp.lfp);
        hold on
        plot(obj.lfp.swr.swr_peaks,ones(length(obj.lfp.swr.swr_on),1),'k*');
        ylabel('broad band')
        ax2=subplot(k,1,2);
        plot(obj.lfp.ts,obj.lfp.swr.swr);
        hold on
        plot(obj.lfp.swr.swr_peaks,ones(length(obj.lfp.swr.swr_on),1),'k*');
        ylabel('SWR')
        ax3=subplot(k,1,3);
        plot(obj.lfp.ts,obj.lfp.delta);
        ylabel('\delta')
        ax4=subplot(k,1,4);
        plot(obj.lfp.ts,obj.lfp.theta);
        ylabel('\theta')
        ax5=subplot(k,1,5);
        plot(obj.lfp.ts,obj.lfp.gamma);
        ylabel('\gamma')
        %         if ~isempty(obj.camera.cam.mvt)
        %             ax6=subplot(k,1,6);
        %             plot(obj.camera.ts_cam,obj.camera.cam.traces);
        %             hold(ax6,'on');
        %             plot(ax6,obj.camera.ts_cam(obj.camera.cam.mvt),obj.camera.cam.traces(obj.camera.cam.mvt),'.');
        %             ylabel('movement')
        %         end
        xlabel('time (sec)')
        linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
        
    case 'deconv_spec'
        if isempty(obj.spec.spectrum)
            error('A spectrum needs to be generated before you can plot it');
        end
        if isempty(obj.twop.deconv)
            error('Deconv must be loaded into current LFP object');
        end
        figure;
        ax1=subplot(2,1,1);
        plot_spec;
        ax2=subplot(2,1,2);
        imagesc(ax2,'xdata',obj.twop.ts,'cdata',fast_smooth(zscore(obj.twop.deconv),2)');
        colormap(ax2,get_colour('black'));
        colorbar;
        linkaxes([ax1 ax2],'x');
        
    case 'channels' %plot all channels extracted from abf, useful for visualizing mixups
        figure;
        k=size(obj.abf.raw,2);
        for i=1:k
            ax(i)=subplot(k,1,i);
            plot(obj.abf.raw(1:1e6,i));
            if ismember(i, obj.abf.Channels)
                ylabel([num2str(i) ' = ' obj.get_channel(i)]);
            else
                ylabel(num2str(i));
            end
        end
        linkaxes(ax,'x');
        
    case 'topography'
        colours=jet;
        bins = size(obj.analysis.stack, 1);
        colours = colours( round(linspace(1, size(colours,1), bins)) , :);
        figure;
        bkg = nan(size(obj.topo.mimg));
        h = imagesc(bkg);
        set(h, 'alphadata', 0);
        colormap jet;
        set(gca, 'xtick',[])
        set(gca, 'ytick',[])
        hold on
        for i=1:bins
            [x,y]=find(obj.topo.loc==i);
            plot(x,y,'.', 'color', colours(i,:));
        end
        axis square
        p = get(gca,'position');
        c = colorbar;
        c.Label.String = 'Position (cm)';
        caxis([0 obj.analysis.vr_length]);
        set(gca, 'position',p);
        
        centroid = false(size(obj.topo.mimg));
        idx = round( obj.topo.centroid ./ obj.topo.FOV .* size(obj.topo.mimg)');
        idx = sub2ind( size(centroid), idx(2,:), idx(1,:) );
        centroid(idx) = true;
        centroid = centroid .* obj.topo.loc;
        centroid(~centroid) = nan;
%         centroid = fast_smooth2(centroid, 'sig', 10);
        centroid = fast_smooth2(centroid, 'method', 'mean', 'wdw', 50);
        centroid = (centroid - 1) ./ (bins - 1) .* obj.analysis.vr_length;
        figure;
        h = imagesc(centroid');
        set(h, 'alphadata', ~isnan(centroid'));
        colormap jet
        axis square
        p = get(gca,'position');
        c = colorbar;
        c.Label.String = 'Position (cm)';
        caxis([0 obj.analysis.vr_length]);
        set(gca, 'position',p);
        
        
end


    function plot_spec()
        if exist('ax1','var')
            imagesc(ax1,'xdata',obj.lfp.spec.t,'ydata',obj.lfp.spec.f,'cdata',obj.lfp.spec.spectrum);
        else
            imagesc('xdata',obj.lfp.spec.t,'ydata',obj.lfp.spec.f,'cdata',obj.lfp.spec.spectrum);
        end
        colormap jet
        xlabel('time (sec)');
        ylabel('frequency (Hz)');
        xlim([obj.lfp.spec.t(1) obj.lfp.spec.t(end)]);
        ylim([obj.lfp.spec.f(1) obj.lfp.spec.f(end)]);
        c=colorbar;
        c.Label.String='log-normalized power';
    end

end