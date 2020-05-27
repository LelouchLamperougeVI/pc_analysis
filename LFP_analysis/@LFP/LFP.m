classdef LFP < handle
    % All the tools for analyzing simultaneous 2p + LFP in one place
    
    properties (GetAccess = 'public', SetAccess = 'protected')
        lfp = struct('lfp',[], 'fs',[], 'ts',[], 'delta',[], 'theta',[], ...
            'gamma',[], 'swr',...
            struct('swr_env',[], 'swr_peaks',[], 'swr_dur',[], ...
            'swr_on',[], 'swr_cyc',[], 'swr',[]) );
        twop = struct( 'fs',[], 'ts',[], 'deconv',[], 'plane', 1, 'numplanes', 1 );
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
%         Channels = [1 2 3 5 5 5]'; %old behavior for RSC RRR (phil trans B)
%         Channels = [1 2 3 6 nan nan]'; %old behavior for RSC RRR (phil trans B)
%         Channels = [1 2 3 6 nan nan 5]'; % licks
        Channels = [1 2 3 5 nan nan nan]'; % Dun topography collab
%         Channels = [1 2 3 5 nan nan]'; %old behavior for RSC RRR (phil trans B)
        %         Channels = [1 2 3 5 5 6]'; %old old behavior for RSC RRR, also works for Ingrid's RRR
        %         Channels = [1 2 3 7 5 5]';
        %         Channels = [1 2 3 5 5 6]';
        %         Channels = [1 2 3 6 8 7]'; %new behavior for RSC RRR
        %         Channels = [1 2 3 4 8 7]'; %new vr behavior for RSC RRR
        %         Channels = [1 3 2 6 5 5]'; % lesion across days
%                 Channels = [1 2 3 4 5 5]'; % aubrey's data
%                 Channels = [1 4 5 6 nan nan]'; % aubrey's data
%         Channels = [1 2 3 5 5 5]'; % haoran's data
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        ChannelNames = {'2p';'chA';'chB';'rwd';'cam';'lfp';'lck'};
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
            [fn,path,mode] = obj.parse_inputs(varargin);
            if isempty(fn)
                [fn,path,mode] = uigetfile({'*.abf','abf file (*.abf)';'behavior.mat','behaviour file (behavior.mat)'});
                switch mode
                    case 1
                        mode = 'abf';
                    case 2
                        mode = 'beh';
                    otherwise
                        error('need to select a valid file');
                end
            end
            
            path = strsplit(path, {'\\','/'});
            if strcmp(mode, 'abf')
                obj.session.date = datetime( regexp(fn, '^\d{4}_(\d{2}_){2}', 'match','once'), 'inputformat', 'yyyy_MM_dd_' );
                obj.session.id = regexp(fn, '_(\d+)\.abf$', 'tokens','once'); obj.session.id = obj.session.id{1};
            else
                obj.session.date = datetime(path, 'inputformat', 'yyy_MM_dd');
                idx = find(~isnat(obj.session.date), 1, 'last');
                obj.session.date = obj.session.date(idx);
                obj.session.id = path{idx+1};
            end
            idx = strcmp( path, datestr(obj.session.date, 'yyyy_mm_dd') );
            path = path(1: find(idx));
            obj.session.animal = path{find(idx)-1};
            if isunix
                path = strjoin( path, '/');
            elseif ispc
                path = strjoin( path, '\');
            end
            obj.session.wd = path;
            
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
            if strcmp(mode, 'abf')
                obj.load(fullfile(path,fn));
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
            if strcmp(mode, 'abf')
                obj.down_sample(obj.default_ops.down_fs);
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
        
        function load(obj, fn) % replace current data with another abf file
            [d,s] = abfload(fn);
            obj.si=s;
            obj.lfp.fs = 1/s * 1e+6;
            obj.raw = d;
            try
                obj.lfp.lfp = d(:,obj.get_channel('lfp'));
            catch
                warning('LFP channel is undefined or unavailable');
            end
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
        
        function perform_analysis(obj,varargin) %Run pc_batch_analysis
%             if ~isempty(obj.analysis)
%                 return
%             end
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
        
        function remove_mvt(obj, mode) % remove moving epochs from deconv detected by camera and belt encoder (whatever is available)
            if nargin < 2
                mode = 'exclude';
            end
            switch mode
                case 'exclude'
                    mode = 0;
                case 'include'
                    mode = 1;
                otherwise
                    error('movement filtration mode unrecognised');
            end
            
            if isempty(obj.twop.deconv)
                error('deconv needs to be loaded into current LFP object first');
            end
            
            if isempty(obj.cam)
                warning('no irCam object loaded; skipping cam movement removal');
            elseif isempty(obj.cam.mvt)
                warning('movement trace hasn''t been extracted from currently loaded irCam object; skipping cam movement removal');
            else
                if ~mode
                    heads=get_head(obj.cam.mvt');
                else
                    heads=get_head(~obj.cam.mvt');
                end
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
                if ~mode
                    thres=abs(obj.behavior.speed_raw)>thres;
                else
                    thres=abs(obj.behavior.speed_raw)<thres;
                end
                obj.twop.deconv(thres,:)=nan;
            else
                thres=noRun(obj.behavior.unit_vel);
                if ~mode
                    thres=abs(obj.behavior.unit_vel)>thres;
                else
                    thres=abs(obj.behavior.unit_vel)<thres;
                end
                obj.twop.deconv(thres,:)=nan;
            end
        end
        
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
        
        function [fn, path, mode] = parse_inputs(obj,inputs)
            mode = 'abf';
            fn = []; path = [];
            count = 1;
            while count <= length(inputs)
                switch(lower(inputs{count}))
                    case {'ch', 'chan', 'channel', 'channels'}
                        obj.set_channels(inputs{count+1});
                    case 'plane'
                        obj.twop.plane = inputs{count+1};
                    otherwise
                        exp = regexp(inputs{count}, {'\.abf$', 'behavior\.mat$'}, 'once');
                        idx = ~cellfun(@isempty, exp);
                        if any(idx)
                            path = strsplit(inputs{count}, {'\\','/'});
                            fn = path{end};
                            if isunix
                                path = strjoin( path(1:end-1), '/');
                            elseif ispc
                                path = strjoin( path(1:end-1), '\');
                            end
                            count = count - 1;
                            
                            switch find(idx)
                                case 1
                                    mode = 'abf';
                                case 2
                                    mode = 'beh';
                            end
                        else
                            error(['argument ''' inputs{count} ''' is undefined']);
                        end
                end
                count = count + 2;
            end
        end
        
        two_photon_ts(obj, deconv)
        irCam_ts(obj)
    end
    
end
