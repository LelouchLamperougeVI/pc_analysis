classdef LFP < handle
    % All the tools for analyzing simultaneous 2p + LFP in one place
    
    properties (GetAccess = 'public', SetAccess = 'protected')
        lfp = struct('lfp',[], 'fs',[], 'ts',[], 'delta',[], 'theta',[], ...
            'gamma',[], 'swr',...
            struct('swr_env',[], 'swr_peaks',[], 'swr_dur',[], ...
            'swr_on',[], 'swr_cyc',[], 'swr',[]), ...
            'spec', struct('spectrum',[],'t',[],'f',[]));
        twop = struct( 'fs',[], 'ts',[], 'deconv',[], 'plane', [], 'numplanes', [] );
        behavior
        camera = struct('ts_cam', [], 'cam', []); % clampex camera pulses and irCam object
        analysis
        abf
        session % session information
        topo = struct('maskNeurons',[], 'mimg',[]);
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        intern
    end
    
    
    methods
        function obj = LFP(varargin) %construct with passed abf file
            if ~isempty(varargin) && iscell(varargin{1})
                varargin = varargin{1};
            end
            obj.set([]);
            obj.set(varargin);
            
            if isempty(obj.intern.fn)
                [obj.intern.fn, obj.intern.path, obj.intern.op_mode] = uigetfile({'*.abf','abf file (*.abf)';'behavior.mat','behaviour file (behavior.mat)'});
                switch obj.intern.op_mode
                    case 1
                        obj.intern.op_mode = 'abf';
                    case 2
                        obj.intern.op_mode = 'beh';
                    otherwise
                        error('need to select a valid file');
                end
            end
            
            obj.intern.path = strsplit(obj.intern.path, {'\\','/'});
            if strcmp(obj.intern.op_mode, 'abf')
                obj.session.date = datetime( regexp(obj.intern.fn, '^\d{4}_(\d{2}_){2}', 'match','once'), 'inputformat', 'yyyy_MM_dd_' );
                obj.session.id = regexp(obj.intern.fn, '_(\d+)\.abf$', 'tokens','once'); obj.session.id = obj.session.id{1};
            else
                obj.session.date = datetime(obj.intern.path, 'inputformat', 'yyy_MM_dd');
                idx = find(~isnat(obj.session.date), 1, 'last');
                obj.session.date = obj.session.date(idx);
                obj.session.id = obj.intern.path{idx+1};
            end
            idx = strcmp( obj.intern.path, datestr(obj.session.date, 'yyyy_mm_dd') );
            obj.intern.path = obj.intern.path(1: find(idx));
            obj.session.animal = obj.intern.path{find(idx)-1};
            if isunix
                obj.intern.path = strjoin( obj.intern.path, '/');
            elseif ispc
                obj.intern.path = strjoin( obj.intern.path, '\');
            end
            obj.session.wd = obj.intern.path;
            
            if exist(fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']), 'file')
                lfp = load(fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']));
                obj = lfp.lfp;
                disp(['loaded existing: ' fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat'])]);
                return
            end
            
            plane_case = dir(fullfile(obj.session.wd, obj.session.id));
            plane_case = regexp( {plane_case([ plane_case.isdir ]).name}, '[Pp]lane\d+', 'match', 'once' );
            plane_case( cellfun(@isempty, plane_case) ) = [];
            plane_case = strrank(plane_case);
            if length(plane_case) > 1
                disp('Multiplane recording detected!!!')
                disp( ['The following planes are available: ' strjoin(plane_case, ' ') ] );
                obj.twop.numplanes = length(plane_case);
            else
                obj.twop.numplanes = 1;
            end
            if ~isempty(plane_case)
                disp( ['loading ' plane_case{obj.twop.plane} '...'] );
                plane = plane_case{obj.twop.plane};
            else
                disp('No planes folder. Loading from session root directory.');
                plane = [];
            end
            
            deconv=load(fullfile(obj.session.wd, obj.session.id, plane, 'deconv.mat'));
            obj.update_channels();
            if strcmp(obj.intern.op_mode, 'abf')
                obj.load_abf(fullfile(obj.intern.path, obj.intern.fn));
                obj.two_photon_ts(deconv.deconv);
            else
                tcs=load(fullfile(obj.session.wd, obj.session.id, plane, 'timecourses.mat'));
                obj.twop.ts = tcs.tcs.tt';
                obj.twop.fs = 1/median(diff(obj.twop.ts));
                try
                    behavior=load(fullfile(obj.session.wd, obj.session.id, plane, 'behavior.mat'));
                catch
                    behavior=load(fullfile(obj.session.wd, obj.session.id, 'behavior.mat'));
                end
                obj.behavior = behavior.behavior;
            end
            try obj.irCam_ts; catch; end
            if strcmp(obj.intern.op_mode, 'abf')
                obj.down_sample(obj.lfp.ops.down_fs);
                obj.filter_bands;
                obj.extract_behaviour;
            end
            
            try
                obj.import_deconv(deconv.deconv);
                disp('Deconv loaded')
            catch
                disp('Failed to autoload deconv :(');
            end
            
            if size(obj.twop.deconv,1) < length(obj.behavior.speed_raw) && ~isempty(obj.twop.deconv)
                disp('resampling behaviour for multiplane...');
                obj.behavior.speed_raw = obj.behavior.speed_raw(obj.twop.plane:obj.twop.numplanes:end);
                if isfield(obj.behavior, 'speed_raw_noSmooth')
                    obj.behavior.speed_raw_noSmooth = obj.behavior.speed_raw_noSmooth(obj.twop.plane:obj.twop.numplanes:end);
                end
            end
            
            if exist(fullfile(obj.session.wd, ['analysis_' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']), 'file')
                analysis=load(fullfile(obj.session.wd, ['analysis_' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']));
                obj.import_analysis(analysis.analysis);
                disp('Analysis loaded');
            elseif exist(fullfile(obj.session.wd, 'analysis.mat'), 'file')
                analysis=load(fullfile(obj.session.wd, 'analysis.mat'));
                obj.import_analysis(analysis.analysis);
                disp('Analysis loaded');
            else
                disp('No analysis file found');
            end
            if exist(fullfile(obj.session.wd, obj.session.id, plane, 'masks_neurons.mat'), 'file')
                maskNeurons = load(fullfile(obj.session.wd, obj.session.id, plane, 'masks_neurons.mat'));
                mimg = load(fullfile(obj.session.wd, obj.session.id, plane, 'mean_img.mat'));
                obj.topo.maskNeurons = maskNeurons.maskNeurons;
                obj.topo.mimg = double(mimg.mimg);
                disp('Topography loaded');
            else
                disp('Failed to load topography');
            end
        end
        
        function import_cam(obj,cam) % import a camera object
            if ~isa(cam,'irCam')
                error('object must be of class ''irCam''');
            end
            obj.camera.cam = cam;
        end
        
        function import_analysis(obj,analysis) % import pc_analysis obtained from run trials
            obj.analysis=analysis;
            obj.analysis.order=get_order(analysis);
        end
        
        function perform_analysis(obj,varargin) %Run pc_batch_analysis
            if isempty(obj.behavior)
                error('no behavioural data currently exists');
            end
            if isempty(obj.twop.deconv)
                error('no deconv data currently exists');
            end
            tcs.tt=obj.twop.ts';
            [beh,dec]=convert_behavior(obj.behavior,tcs,obj.twop.deconv);
            obj.analysis=pc_batch_analysis(beh,dec,varargin);
            obj.analysis.order=get_order(obj.analysis);
        end
        
        function chan = get_channel(obj,ch) % return currently listed channels
            if isa(ch,'char')
                chan = obj.abf.chan_tb(ch,:).Variables;
            elseif isa(ch,'double')
                chan=obj.abf.chan_tb(obj.abf.chan_tb.Variables==ch,:).Row{1};
                %                 chan = [num2str(obj.chan(ch,:).Variables) ' = ' obj.chan(ch,:).Row{1}];
            else
                error('unrecognized input');
            end
        end
        
        function reset(obj) % reset to original lfp from abf
            obj.lfp.fs = 1/obj.abf.si * 1e+6;
            obj.lfp.lfp = obj.abf.raw(:,obj.get_channel('lfp'));
        end
        
        function invert_pol(obj) % invert the polarities between electrode tips
            m=median(obj.lfp.lfp);
            obj.lfp.lfp=-(obj.lfp.lfp-m)+m;
        end
        
        remove_mvt(obj, mode);
        import_deconv(obj,deconv);
        down_sample(obj,Fs);
        spectrum(obj,win,range);
        detrend(obj,wdw);
        reference60(obj,len);
        filter_bands(obj);
        extract_behaviour(obj);
        enregistrer(obj);
        topography(obj, FOV);
    end
    
    methods (Access = private)
        function update_channels(obj)
            if min(size(obj.abf.Channels)) > 1
                error('Channels must be a 1D vector');
            end
            if length(obj.abf.Channels) < length(obj.abf.ChannelNames) || length(obj.abf.Channels) < length(obj.abf.ChannelNames)
                error(['Channels must be defined as a vector of ' num2str(length(obj.abf.ChannelNames)) ' integers']);
            end
            if size(obj.abf.Channels,2)>1
                obj.abf.Channels = obj.abf.Channels';
            end
            obj.abf.chan_tb = table(obj.abf.Channels,'RowNames',obj.abf.ChannelNames);
        end
        
        function load(obj, fn)
            [d,s] = abfload(fn);
            obj.abf.si=s;
            obj.lfp.fs = 1/s * 1e+6;
            obj.abf.raw = d;
            try
                obj.lfp.lfp = d(:,obj.get_channel('lfp'));
            catch
                warning('LFP channel is undefined or unavailable');
            end
        end
        
        two_photon_ts(obj, deconv)
        irCam_ts(obj)
    end
    
end
