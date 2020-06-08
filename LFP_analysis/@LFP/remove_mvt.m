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

if isempty(obj.camera.cam)
    warning('no irCam object loaded; skipping cam movement removal');
elseif isempty(obj.camera.cam.mvt)
    warning('movement trace hasn''t been extracted from currently loaded irCam object; skipping cam movement removal');
else
    if ~mode
        heads=get_head(obj.camera.cam.mvt');
    else
        heads=get_head(~obj.camera.cam.mvt');
    end
    heads=obj.camera.ts_cam(heads);
    tails=get_head(obj.camera.cam.mvt(end:-1:1)');
    tails=obj.camera.ts_cam(tails(end:-1:1));
    
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