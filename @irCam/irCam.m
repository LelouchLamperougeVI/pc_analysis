classdef irCam < handle
    % workflow: add_mask (optional - add_label) --> extract_traces --> filter_traces --> detect_mvt --> player --> export_mvt
    properties
        frames
        traces
        mvt
        dims=zeros(1,2);
        fs
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        cam
        original_traces
        trace_std
        num_frames
        nose_mask
        limb_mask
        
        masks
        mask_labels = {'nose',...
                       'limb' };
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        default_labels_idx
    end
    
    methods
        function obj = irCam(fn)
            if nargin < 1
                [fn, path] = uigetfile('*.avi');
                fn = fullfile(path,fn);
            end
            mat_fn = [fn(1:end-3) 'mat'];
            if exist(mat_fn, 'file')
                load(mat_fn);
                
            else
                obj.cam=VideoReader(fn);
%                 save(mat_fn, 'obj.cam');
            end
            obj.dims(1)=obj.cam.Height;
            obj.dims(2)=obj.cam.Width;
            obj.num_frames=obj.cam.FrameRate*obj.cam.Duration;
            obj.fs = obj.num_frames / obj.cam.Duration;
            
            obj.default_labels_idx = length(obj.mask_labels);
            obj.masks = false(obj.dims(1), obj.dims(2), obj.default_labels_idx);
            
%             obj.choose_nose;
%             obj.extract_traces;
%             obj.remove_steps;
%             obj.detect_mvt;
        end
        
        function add_mask(obj,label)
            % choose an ROI for a label
            mask = obj.gui_get_bounds;
            close gcf
            if nargin < 2
                label = listdlg('liststring',obj.mask_labels, 'selectionmode','single', 'promptstring','Select mask label:');
                if isempty(label)
                    error('You must choose a label! If you need a custom label, you may add one with ''add_label''');
                end
            else
                label = cellfun(@(x) strcmpi(label, x), obj.mask_labels);
                if ~sum(label)
                    error('The label name you have passed to this function is not recognized. If you need a custom label, you may add one with ''add_label''');
                end
                label = find(label);
            end
            
            obj.masks(:,:,label) = mask;
        end
        
        function add_label(obj, label)
            % add a new label
            if nargin < 2
                error('You must provide a name for you custom label')
            end
            if ~isa(label, 'char')
                error('Your label must be provided as a string')
            end
            
            obj.mask_labels{end+1} = label;
            obj.masks(:,:,end+1) = false(size(obj.masks,1),size(obj.masks,2));
        end
        
        function remove_label(obj, label)
            % remove a label from the list
            idx = cellfun(@(x) strcmpi(label, x), obj.mask_labels);
            if ~sum(idx)
                error('Dummy, the label you''re trying to remove doesn''t appear to exist in the first place...')
            end
            if idx <= obj.default_labels_idx
                error('Sorry, I can''t let you remove a label that is part of the default list of masks. You can use clear_mask instead so that its trace wouldn''t be processed.')
            end
            obj.mask_labels(idx) = [];
            obj.masks(:,:,idx) = [];
        end
        
        function clear_mask(obj, mask)
            % clear the ROI selection for a label
            switch class(mask)
                case {'double', 'logical'}
                    obj.get_label(mask); %check for error
                case 'char'
                    mask = get_label(mask);
                otherwise
                    error('Invalid variable type passed as input')
            end
            obj.masks(:,:, mask) = false;
        end
        
        function rval = get_label(obj, label)
            %Get the name or index of a mask label
            switch class(label)
                case {'double','logical'}
                    if (length(label) > length(obj.mask_labels)) || any(label > length(obj.mask_labels))
                       error('The input index exceeds the range of available labels') 
                    end
                    rval = strjoin(obj.mask_labels(label), ', ');
                case 'char'
                    idx = cellfun(@(x) strcmpi(label, x), obj.mask_labels);
                    if ~sum(idx)
                        error('The search term was not found to be part of the labels list')
                    end
                    rval = find(idx);
                otherwise
                    error('The input is fucked...')
            end
        end
        
        function choose_nose(obj)
            % DEPRECATED
            obj.nose_mask=obj.gui_get_bounds;
        end
        
        function choose_limb(obj)
            % DEPRECATED
            obj.limb_mask=obj.gui_get_bounds;
        end
        
        function detect_mvt(obj,sig,dilation,gap,dur)
            % detect motion artifacts
            % pipes: MAD --> dilate --> fill_gaps --> duration filt
            % sig: number of MADs from median (default 3)
            % dilation: dilate the detected regions (default .5 sec - two-tailed; i.e. not 1 sec each tails)
            % gap: number of samples between gaps to fill in (default 2*dilation)
            % dur: minimum duration of movement epoch (default 2 sec)
            if nargin<2; sig=3; end
            if nargin<3; dilation=.5; end
            if nargin<4; gap=2*dilation; end
            if nargin<5; dur=2; end %dilate the detection
            gap=obj.cam.FrameRate*gap;
            obj.mvt = any( abs(obj.traces)> (sig .* mad(obj.traces)) ,2);
            obj.mvt = logical( conv2(obj.mvt, ones(obj.cam.FrameRate*dilation,1), 'same') );
            obj.mvt = fill_gaps(obj.mvt,gap);
            obj.mvt = double(obj.mvt);
            
            heads = find(get_head(obj.mvt));
            tails = sort( length(obj.mvt)+1 - find( get_head( obj.mvt(end:-1:1) ) ) );
            idx = find((tails - heads +1) < (dur*obj.cam.FrameRate));
            for i=1:length(idx)
                obj.mvt(heads(idx(i)) : tails(idx(i))) = 2;
            end
            
        end
        
        function mvt = export_mvt(obj)
            % export a logical array for timestamps where movements occured
            mvt = obj.mvt;
            mvt(mvt==2) = 0;
            mvt = logical(mvt);
        end
        
        extract_traces(obj);
        remove_steps(obj,thres,samples);
        filter_traces(obj,freq);
        player(obj);
    end
    
    methods(Access = private)
        mask = gui_get_bounds(obj);
    end
end