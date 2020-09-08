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
    thres = mk_sync(thres, length(obj.twop.planes.planes));
    obj.twop.deconv(thres,:)=nan;
else
    thres=noRun(obj.behavior.unit_vel);
    if ~mode
        thres=abs(obj.behavior.unit_vel)>thres;
    else
        thres=abs(obj.behavior.unit_vel)<thres;
    end
    thres = mk_sync(thres, length(obj.twop.planes.planes));
    obj.twop.deconv(thres,:)=nan;
end

obj.filter_bands;


function mvt = mk_sync(mvt, numPlanes)
% for multi-plane data, make sure that the number of frames to be removed
% is a multiple of the number of planes analysed

l = find(get_head(mvt));
r = get_head(mvt(end:-1:1));
r = find(r(end:-1:1));
residual = mod(r - l + 1, numPlanes);
residual(~~residual) = numPlanes - residual(~~residual);
% MATLAB, for the love of fuck, please make it possible to have do-while loops...
while any(~~residual)
    for ii = find(~~residual)'
        lidx = l(ii) - ceil(residual(ii) / 2);
        ridx = r(ii) + floor(residual(ii) / 2);
        if lidx < 1
            mvt(1 : l(ii)) = 1;
            mvt(r(ii) : ridx + (1 - lidx)) = 1;
            continue
        end
        
        if ridx > length(mvt)
            mvt(r(ii) : end) = 1;
            mvt(lidx - (ridx - length(mvt)) : l(ii)) = 1;
            continue
        end
        
        mvt(lidx : l(ii)) = 1;
        mvt(r(ii) : ridx) = 1;
    end
    
    l = find(get_head(mvt));
    r = get_head(mvt(end:-1:1));
    r = find(r(end:-1:1));
    residual = mod(r - l + 1, numPlanes);
    residual(~~residual) = numPlanes - residual(~~residual);
end