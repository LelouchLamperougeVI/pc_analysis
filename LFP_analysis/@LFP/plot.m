function plot(obj,option)
% Options:
%   'spectrum'
%   'bands'
%   'deconv_spec'

switch option
    case 'spectrum'
        if isempty(obj.spec.spectrum)
            error('A spectrum needs to be generated before you can plot it');
        end
        figure;
        plot_spec;
        
    case 'bands'
        k=5;
        %         if ~isempty(obj.cam.mvt)
        %             k=k+1;
        %         end
        t=linspace(0,length(obj.lfp)/obj.fs,length(obj.lfp));
        figure;
        ax1=subplot(k,1,1);
        plot(t,obj.lfp);
        ylabel('broad band')
        ax2=subplot(k,1,2);
        plot(t,obj.delta);
        ylabel('\delta')
        ax3=subplot(k,1,3);
        plot(t,obj.theta);
        ylabel('\theta')
        ax4=subplot(k,1,4);
        plot(t,obj.gamma);
        ylabel('\gamma')
        ax5=subplot(k,1,5);
        plot(t,obj.swr);
        ylabel('SWR')
        %         if ~isempty(obj.cam.mvt)
        %             ax6=subplot(k,1,6);
        %             plot(obj.ts_cam,obj.cam.traces);
        %             hold(ax6,'on');
        %             plot(ax6,obj.ts_cam(obj.cam.mvt),obj.cam.traces(obj.cam.mvt),'.');
        %             ylabel('movement')
        %         end
        xlabel('time (sec)')
        linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
        
    case 'deconv_spec'
        if isempty(obj.spec.spectrum)
            error('A spectrum needs to be generated before you can plot it');
        end
        if isempty(obj.deconv)
            error('Deconv must be loaded into current LFP object');
        end
        figure;
        ax1=subplot(2,1,1);
        plot_spec;
        ax2=subplot(2,1,2);
        imagesc(ax2,'xdata',obj.ts_2p,'cdata',fast_smooth(zscore(obj.deconv),2)');
        colormap(ax2,get_colour('black'));
        colorbar;
        linkaxes([ax1 ax2],'x');
        
end


    function plot_spec()
        if exist('ax1','var')
            imagesc(ax1,'xdata',obj.spec.t,'ydata',obj.spec.f,'cdata',obj.spec.spectrum);
        else
            imagesc('xdata',obj.spec.t,'ydata',obj.spec.f,'cdata',obj.spec.spectrum);
        end
        colormap jet
        xlabel('time (sec)');
        ylabel('frequency (Hz)');
        xlim([obj.spec.t(1) obj.spec.t(end)]);
        ylim([obj.spec.f(1) obj.spec.f(end)]);
        c=colorbar;
        c.Label.String='log-normalized power';
    end

end