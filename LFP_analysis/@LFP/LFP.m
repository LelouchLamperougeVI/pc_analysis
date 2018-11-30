classdef LFP < handle
    properties (GetAccess = 'public', SetAccess = 'private')
        lfp %unfiltered lfp trace
        fs %sampling rate in Hz
        fs_2p
        ts_2p
        ts_cam
        cam
        delta
        theta
        gamma
        swr
        spec = struct('spectrum',[],'t',[],'f',[]);
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        ChannelNames = {'2p';'chA';'chB';'cam';'lfp'};
        Channels = [1 2 3 8 7]';
        chan;
        raw
        f
        f60_env = 1;
    end
    
    properties (GetAccess = 'private', Constant)
        nfft = 2^16;
    end
    
    methods
        function obj = LFP(fn) %construct with passed abf file
            obj.update_channels();
            if nargin > 0
                obj.load(fn);
                obj.two_photon_ts;
                try
                    obj.irCam_ts;
                catch
                end
            end
        end
        
        function load(obj, fn)
            [d,si] = abfload(fn);
            obj.fs = 1/si * 1e+6;
            obj.raw = d;
            obj.lfp = d(:,obj.get_channel('lfp'));
        end
        
        function import_cam(obj,cam)
            if ~isa(cam,'irCam')
                error('object must be of class ''irCam''');
            end
            obj.cam=cam;
        end
        
        function set_channels(obj, ch)
            obj.Channels = ch;
            obj.update_channels();
        end
        
        function chan = get_channel(obj,ch)
            chan = obj.chan(ch,:).Variables;
        end
        
        function update(obj)
        end
        
        down_sample(obj,Fs);
        spectrum(obj,win);
        detrend(obj,wdw);
        reference60(obj,len);
        filter_bands(obj);
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