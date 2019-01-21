classdef LFP < handle
    properties (GetAccess = 'public', SetAccess = 'private')
        lfp %unfiltered lfp trace
        fs %sampling rate in Hz
        t %timestamps series
        fs_2p
        ts_2p
        behavior
        ts_cam
        cam % irCam object
        deconv
        analysis
        ensemble % basic_ensemble object
        delta
        theta
        gamma
        swr
        swr_env % ripple power/amplitude envelope
        swr_peaks %ripple peak onset times
        swr_dur %ripple durations
        swr_on %ripple onset times
        swr_cyc %number of ripple frequency cycles
        
        spec = struct('spectrum',[],'t',[],'f',[]);
        
        %         Channels = [1 2 3 5 8 7]';
        %         Channels = [1 2 3 5 8 6]';
        Channels = [1 2 3 7 8 5]';
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        ChannelNames = {'2p';'chA';'chB';'rwd';'cam';'lfp'};
        chan;
        raw
        si
        f
        f60_env = 1;
    end
    
    properties (GetAccess = 'private', Constant)
        nfft = 2^16;
        default_ops = struct('down_fs',2e3,'spec_wdw',.5);
    end
    
    methods
        function obj = LFP(fn) %construct with passed abf file
            obj.update_channels();
            obj.load(fn);
            obj.two_photon_ts;
            try
                obj.irCam_ts;
            catch
            end
            obj.down_sample(obj.default_ops.down_fs);
            obj.filter_bands;
            obj.extract_behaviour;
        end
        
        function load(obj, fn) % replace current data with another abf file
            [d,s] = abfload(fn);
            obj.si=s;
            obj.fs = 1/s * 1e+6;
            obj.raw = d;
            obj.lfp = d(:,obj.get_channel('lfp'));
        end
        
        function import_cam(obj,cam) % import a camera object
            if ~isa(cam,'irCam')
                error('object must be of class ''irCam''');
            end
            obj.cam=cam;
        end
        
        function import_deconv(obj,deconv) % import deconv
            if size(deconv,1)>18000
                obj.deconv=stupid_windows_fs(deconv);
            else
                obj.deconv=deconv;
            end
        end
        
        function import_analysis(obj,analysis) % import pc_analysis obtained from run trials
            obj.analysis=analysis;
            obj.analysis.order=get_order(analysis);
        end
        
        function import_behaviour(obj,behavior)
            obj.behavior=behavior;
        end
        
        function perform_analysis(obj,varargin) %Run pc_batch_analysis
            if isempty(obj.behavior)
                error('no behavioural data currently exists');
            end
            if isempty(obj.deconv)
                error('no deconv data currently exists');
            end
            tcs.tt=obj.ts_2p;
            [beh,dec]=convert_behavior(obj.behavior,tcs,obj.deconv);
            obj.analysis=pc_batch_analysis(beh,dec,varargin);
            obj.analysis.order=get_order(obj.analysis);
        end
        
        function detect_sce(obj,varargin) % construct basic_ensemble object 
            if isempty(obj.deconv)
                error('no deconv data available');
            end
            obj.ensemble = basic_ensemble(obj.deconv,obj.ts_2p,varargin);
            if ~isempty(obj.analysis)
                obj.ensemble.order=get_order(obj.analysis);
            end
            obj.ensemble.lfp=obj;
        end
        
        function set_channels(obj, ch) % change channels configuration from default
            obj.Channels = ch;
            obj.update_channels();
        end
        
        function chan = get_channel(obj,ch) % return currently listed channels
            if isa(ch,'char')
                chan = obj.chan(ch,:).Variables;
            elseif isa(ch,'double')
                chan=obj.chan(obj.chan.Variables==ch,:).Row{1};
%                 chan = [num2str(obj.chan(ch,:).Variables) ' = ' obj.chan(ch,:).Row{1}];
            else
                error('unrecognized input');
            end
        end
        
        function reset(obj) % reset to original lfp from abf
            obj.fs = 1/obj.si * 1e+6;
            obj.lfp = obj.raw(:,obj.get_channel('lfp'));
        end
        
        function invert_pol(obj) % invert the polarities between electrode tips
            m=median(obj.lfp);
            obj.lfp=-(obj.lfp-m)+m;
        end
        
        function remove_mvt(obj) % remove moving epochs from deconv detected by camera and belt encoder (whatever is available)
            if isempty(obj.deconv)
                error('deconv needs to be loaded into current LFP object first');
            end
            
            if isempty(obj.cam)
                warning('no irCam object loaded; skipping cam movement removal');
            elseif isempty(obj.cam.mvt)
                warning('movement trace hasn''t been extracted from currently loaded irCam object; skipping cam movement removal');
            else
                heads=get_head(obj.cam.mvt');
                heads=obj.ts_cam(heads);
                tails=get_head(obj.cam.mvt(end:-1:1)');
                tails=obj.ts_cam(tails(end:-1:1));
                
                heads(tails<obj.ts_2p(1))=[];tails(tails<obj.ts_2p(1))=[];
                tails(heads>obj.ts_2p(end))=[];heads(heads>obj.ts_2p(end))=[];
                
                for i=1:length(heads)
                    idx=[find(obj.ts_2p>heads(i),1) find(obj.ts_2p>tails(i),1)];
                    idx(1)=~(idx(1)<1)*idx(1) + 1;
                    idx(2)=~(idx(2)>length(obj.ts_2p)).*idx(2) + (idx(2)>length(obj.ts_2p))*length(obj.ts_2p);
                    obj.deconv(idx(1):idx(2),:)=nan;
                end
            end
            
            if isempty(obj.behavior)
                warning('no behavioural data loaded; skipping encoder movement removal');
            elseif ~isfield(obj.behavior,'unit_vel')
                thres=noRun(obj.behavior.speed_raw);
                thres=abs(obj.behavior.speed_raw)>thres;
                obj.deconv(thres,:)=nan;
            else
                thres=noRun(obj.behavior.unit_vel);
                thres=abs(obj.behavior.unit_vel)>thres;
                obj.deconv(thres,:)=nan;
            end
        end
        
        down_sample(obj,Fs);
        spectrum(obj,win,range);
        detrend(obj,wdw);
        reference60(obj,len);
        filter_bands(obj);
        extract_behaviour(obj);
    end
    
    methods (Access = private)
        function update_channels(obj)
            if min(size(obj.Channels)) > 1
                error('Channels must be a 1D vector');
            end
            if length(unique(obj.Channels)) < length(obj.ChannelNames) || length(obj.Channels) < length(obj.ChannelNames)
                error(['Channels must be defined as a vector of ' num2str(length(obj.ChannelNames)) ' unique integers']);
            end
            if size(obj.Channels,2)>1
                obj.Channels = obj.Channels';
            end
            obj.chan = table(obj.Channels,'RowNames',obj.ChannelNames);
        end
        
        two_photon_ts(obj)
        irCam_ts(obj)
    end
    
end