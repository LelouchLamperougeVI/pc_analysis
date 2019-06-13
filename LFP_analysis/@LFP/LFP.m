classdef LFP < handle
% All the tools for analyzing simultaneous 2p + LFP in one place

    properties (GetAccess = 'public', SetAccess = 'protected')
        lfp = struct('lfp',[], 'fs',[], 'ts',[], 'delta',[], 'theta',[], ...
                    'gamma',[], 'swr',...
                    struct('swr_env',[], 'swr_peaks',[], 'swr_dur',[], ...
                    'swr_on',[], 'swr_cyc',[], 'swr',[]) );
        twop = struct( 'fs',[], 'ts',[], 'deconv',[] );
        behavior
        ts_cam
        cam % irCam object
        analysis
        ensemble % basic_ensemble object
        
        lfp_mvt %500-1k Hz band signal for movement detection
        
        spec = struct('spectrum',[],'t',[],'f',[]);
        
        session % session information
        topo = struct('maskNeurons',[], 'mimg',[]);
        
%         Channels = [1 2 3 6 8 7]';
%         Channels = [1 2 3 4 5 6]';
%         Channels = [1 2 3 5 4 6]';
%         Channels = [1 2 3 5 8 6]';
%         Channels = [1 2 3 7 8 5]';
%         Channels = [1 2 3 6 5 5]';
        Channels = [1 2 3 5 5 5]'; %old behavior for RSC RRR
%         Channels = [1 2 3 5 5 6]'; %old old behavior for RSC RRR, also works for Ingrid's RRR
%         Channels = [1 2 3 7 5 5]';
%         Channels = [1 2 3 5 5 6]';
%         Channels = [1 2 3 6 8 7]'; %new behavior for RSC RRR
%         Channels = [1 2 3 4 8 7]'; %new vr behavior for RSC RRR
%         Channels = [1 3 2 6 5 5]';
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        ChannelNames = {'2p';'chA';'chB';'rwd';'cam';'lfp'};
        chan;
        raw
        si
        f
        f60_env = 1;
    end
    
    properties (GetAccess = 'protected', Constant)
        nfft = 2^16;
        default_ops = struct('down_fs',2e3,'spec_wdw',.5, 'swr_cyc',3, 'swr_gap', .25, 'FOV', 835.76 .* ones(2,1));
    end
    
    methods
        function obj = LFP(varargin) %construct with passed abf file
            if ~isempty(varargin) && iscell(varargin{1})
                varargin = varargin{1};
            end
            if isempty(varargin)
                [fn,path]=uigetfile('*.abf');
            else
                fn = varargin{1};
            	path = strsplit(fn, {'\\','/'});
                fn = path{end};
                if isunix
                    path = strjoin( path(1:end-1), '/');
                elseif ispc
                    path = strjoin( path(1:end-1), '\');
                end
            end
            obj.session.date = datetime( regexp(fn, '^\d{4}_(\d{2}_){2}', 'match','once'), 'inputformat', 'yyyy_MM_dd_' );
            obj.session.id = regexp(fn, '_(\d+).abf$', 'tokens','once'); obj.session.id = obj.session.id{1};
            path = strsplit(path, {'\\','/'});
            idx = strcmp( path, datestr(obj.session.date, 'yyyy_mm_dd') );
            path = path(1: find(idx));
            obj.session.animal = path{find(idx)-1};
            if isunix
                path = strjoin( path, '/');
            elseif ispc
                path = strjoin( path, '\');
            end
            obj.session.wd = path;
            
            obj.update_channels();
            obj.load(fullfile(path,fn));
            obj.two_photon_ts;
            try obj.irCam_ts; catch; end
            obj.down_sample(obj.default_ops.down_fs);
            obj.filter_bands;
            obj.extract_behaviour;
            
            try
                deconv=load(fullfile(obj.session.wd, obj.session.id, 'plane1', 'deconv.mat'));
                obj.import_deconv(deconv.deconv);
                disp('Deconv loaded')
            catch
                disp('Failed to autoload deconv :(');
            end
            if exist(fullfile(obj.session.wd, 'analysis.mat'), 'file')
                analysis=load(fullfile(obj.session.wd, 'analysis.mat'));
                obj.import_analysis(analysis.analysis);
                disp('Analysis loaded');
            else
                disp('No analysis file found');
            end
            if exist(fullfile(obj.session.wd, obj.session.id, 'plane1', 'masks_neurons.mat'), 'file')
                maskNeurons = load(fullfile(obj.session.wd, obj.session.id, 'plane1', 'masks_neurons.mat'));
                mimg = load(fullfile(obj.session.wd, obj.session.id, 'plane1', 'mean_img.mat'));
                obj.topo.maskNeurons = maskNeurons.maskNeurons;
                obj.topo.mimg = double(mimg.mimg);
                disp('Topography loaded');
            else
                disp('Failed to load topography');
            end
        end
        
        function load(obj, fn) % replace current data with another abf file
            [d,s] = abfload(fn);
            obj.si=s;
            obj.lfp.fs = 1/s * 1e+6;
            obj.raw = d;
            obj.lfp.lfp = d(:,obj.get_channel('lfp'));
        end
        
        function import_cam(obj,cam) % import a camera object
            if ~isa(cam,'irCam')
                error('object must be of class ''irCam''');
            end
            obj.cam=cam;
        end
        
        function import_deconv(obj,deconv) % import deconv
            if exist(fullfile(obj.session.wd, 'fkd18k'), 'file')
                fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'r');
                if fscanf(fid, '%d')
                    obj.twop.deconv = stupid_windows_fs(deconv);
                    disp('18k correction applied');
                else
                    obj.twop.deconv = deconv;
                    if size(deconv,1)>18000
                        disp('18k correction unnecessary');
                    end
                end
                fclose(fid);
                return
            end
            if size(deconv,1)>18000 && ~isempty(obj.behavior)
                tcs.tt=obj.twop.ts';
                [vel,temp] = convert_behavior(obj.behavior,tcs,deconv);
                vel = vel.unit_vel;
                temp = fast_smooth(temp,obj.twop.fs*.5);
                if size(temp,1) < length(vel); vel=vel(1:size(temp,1)); end
                ori_r = corr(temp, vel');
                
                [vel,temp] = convert_behavior(obj.behavior,tcs,stupid_windows_fs(deconv));
                vel = vel.unit_vel;
                temp = fast_smooth(temp,obj.twop.fs*.5);
                if size(temp,1) < length(vel); vel=vel(1:size(temp,1)); end
                alt_r = corr(temp, vel');

                if prctile(abs(alt_r), 70) > prctile(abs(ori_r), 70)
                    disp('DANGER! A rupture in Space-Time continuum has been detected!');
                    uinp = input('Should I fix it for you? Y/N [Y] ', 's');
                    if strcmpi(uinp,'y') || strcmpi(uinp,'yes') || isempty(uinp)
                        obj.twop.deconv=stupid_windows_fs(deconv);
                        disp('OK. I''ve fixed it for you. But know that this only works for behavioural data. It''s up to you to remember to correct it during rest.');                        
                        fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'w');
                        fprintf(fid, '%d', 1);
                        fclose(fid);
                        return
                    end
                    disp('I surely hope you know what you''re doing...');
                end
            end
            obj.twop.deconv=deconv;
            fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'w');
            fprintf(fid, '%d', 0);
            fclose(fid);
        end
        
        function import_analysis(obj,analysis) % import pc_analysis obtained from run trials
            obj.analysis=analysis;
            obj.analysis.order=get_order(analysis);
        end
        
        function import_behaviour(obj,behavior) % import Dun's behavior format (can use extract_behavior instead)
            obj.behavior=behavior;
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
        
%         function detect_sce(obj,varargin) % construct basic_ensemble object 
%             if isempty(obj.twop.deconv)
%                 error('no deconv data available');
%             end
%             obj.ensemble = basic_ensemble(obj.twop.deconv,obj.twop.ts,varargin);
%             if ~isempty(obj.analysis)
%                 obj.ensemble.pc_order=get_order(obj.analysis);
%             end
%             obj.ensemble.lfp=obj;
%         end
        
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
            obj.lfp.fs = 1/obj.si * 1e+6;
            obj.lfp.lfp = obj.raw(:,obj.get_channel('lfp'));
        end
        
        function invert_pol(obj) % invert the polarities between electrode tips
            m=median(obj.lfp.lfp);
            obj.lfp.lfp=-(obj.lfp.lfp-m)+m;
        end
        
        function remove_mvt(obj) % remove moving epochs from deconv detected by camera and belt encoder (whatever is available)
            if isempty(obj.twop.deconv)
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
                
                heads(tails<obj.twop.ts(1))=[];tails(tails<obj.twop.ts(1))=[];
                tails(heads>obj.twop.ts(end))=[];heads(heads>obj.twop.ts(end))=[];
                
                for i=1:length(heads)
                    idx=[find(obj.twop.ts>heads(i),1) find(obj.twop.ts<tails(i),1,'last')];
                    idx(1)=~(idx(1)<1)*idx(1) + (idx(1)<1);
                    idx(2)=~(idx(2)>length(obj.twop.ts)).*idx(2) + (idx(2)>length(obj.twop.ts))*length(obj.twop.ts);
                    obj.twop.deconv(idx(1):idx(2),:)=nan;
                end
            end
            
            if isempty(obj.behavior)
                warning('no behavioural data loaded; skipping encoder movement removal');
            elseif ~isfield(obj.behavior,'unit_vel')
                thres=noRun(obj.behavior.speed_raw);
                thres=abs(obj.behavior.speed_raw)>thres;
                obj.twop.deconv(thres,:)=nan;
            else
                thres=noRun(obj.behavior.unit_vel);
                thres=abs(obj.behavior.unit_vel)>thres;
                obj.twop.deconv(thres,:)=nan;
            end
        end
        
        down_sample(obj,Fs);
        spectrum(obj,win,range);
        detrend(obj,wdw);
        reference60(obj,len);
        filter_bands(obj);
        extract_behaviour(obj);
        enregistrer(obj, fn);
        topography(obj, FOV);
    end
    
    methods (Access = private)
        function update_channels(obj)
            if min(size(obj.Channels)) > 1
                error('Channels must be a 1D vector');
            end
            if length(obj.Channels) < length(obj.ChannelNames) || length(obj.Channels) < length(obj.ChannelNames)
                error(['Channels must be defined as a vector of ' num2str(length(obj.ChannelNames)) ' integers']);
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