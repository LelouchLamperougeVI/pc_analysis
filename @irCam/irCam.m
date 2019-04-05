classdef irCam < handle
    properties
        frames
        traces
        mvt
        dims=zeros(1,2);
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        cam
        original_traces
        num_frames
        nose_mask
        limb_mask
    end
    
    methods
        function obj = irCam(fn)
            if nargin < 1
                [fn, path] = uigetfile('*.avi');
                fn = fullfile(path,fn);
            end
            obj.cam=VideoReader(fn);
            obj.dims(1)=obj.cam.Height;
            obj.dims(2)=obj.cam.Width;
            obj.num_frames=obj.cam.FrameRate*obj.cam.Duration;
            
            obj.choose_nose;
            obj.extract_traces;
            obj.remove_steps;
            obj.detect_mvt;
        end
        
        function choose_nose(obj)
            obj.nose_mask=obj.gui_get_bounds;
        end
        
        function choose_limb(obj)
            obj.limb_mask=obj.gui_get_bounds;
        end
        
        function detect_mvt(obj,sig,gap)
            % detect motion artifacts
            % sig: number of SDs from mean (default 3)
            % gap: number of samples between gaps to fill in (default 3 sec)
            if nargin<2; sig=3; end
            if nargin<3; gap=3; end
            gap=obj.cam.FrameRate*gap;
            obj.mvt=abs(obj.traces-mean(obj.traces))>(std(obj.traces)*sig/2);
            obj.mvt=fill_gaps(obj.mvt,gap);
        end
        
        extract_traces(obj);
        remove_steps(obj,thres,samples)
        player(obj);
    end
    
    methods(Access = private)
        mask = gui_get_bounds(obj);
    end
end