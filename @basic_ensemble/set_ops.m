function set_ops(obj,inputs)
if isempty(obj.ops)
    ops.sig=.05;
    ops.thres=3;
    ops.off_thres=1;
    ops.gaps=.25;
    
    ops.freq=[0 400];
    ops.wdw=[-.5 .5];
    ops.wdw_size=.01;
    
    obj.ops=ops;
end

idx=1;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            obj.ops.sig=inputs{idx};
        case 'thres'
            idx=idx+1;
            obj.ops.thres=inputs{idx};
        case 'gaps'
            idx=idx+1;
            obj.ops.gaps=inputs{idx};
        case 'freq'
            idx=idx+1;
            obj.ops.freq=inputs{idx};
        case 'wdw'
            idx=idx+1;
            obj.ops.wdw=inputs{idx};
        case 'wdw_size'
            idx=idx+1;
            obj.ops.wdw_size=inputs{idx};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx=idx+1;
end