function set_ops(obj,varargin)
if isempty(obj.ops) || nargin<2
    ops.sig=.05;
    ops.thres=3;
    ops.off_thres=1;
    ops.gaps=.25;
    
    ops.freq=[0 400];
    ops.wdw=[-.5 .5];
    ops.wdw_size=.01;
    
    ops.e_size=10;
    ops.e_prctile=.1;
    
    ops.sce_dur=[0 .8]; % 0 to 800 ms SCE events
    
    obj.ops=ops;
    if nargin<2
        return;
    end
end

idx=1;
while(idx<length(varargin))
    switch lower(varargin{idx})
        case 'sig'
            idx=idx+1;
            obj.ops.sig=varargin{idx};
        case 'thres'
            idx=idx+1;
            obj.ops.thres=varargin{idx};
        case 'gaps'
            idx=idx+1;
            obj.ops.gaps=varargin{idx};
            
        case 'freq'
            idx=idx+1;
            obj.ops.freq=varargin{idx};
        case 'wdw'
            idx=idx+1;
            obj.ops.wdw=varargin{idx};
        case 'wdw_size'
            idx=idx+1;
            obj.ops.wdw_size=varargin{idx};
            
        case 'e_size'
            idx=idx+1;
            obj.ops.e_size=varargin{idx};
        case 'e_prctile'
            idx=idx+1;
            obj.ops.e_prctile=varargin{idx};
            
        case 'sce_dur'
            idx=idx+1;
            obj.ops.sce_dur=varargin{idx};
        otherwise
            error(['''' varargin{idx} ''' is not a valid parameter']);
    end
    idx=idx+1;
end