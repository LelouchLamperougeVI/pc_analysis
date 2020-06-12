function two_photon_ts(obj, deconv)
% extract 2p timestamps

ts = obj.abf.raw(:,obj.get_channel('2p')) < 3;
ts = get_head(ts);
ts = (find(ts) - 1)./obj.lfp.fs;
if ts(1) == 0
    ts(1) = [];
end

obj.twop.fs = 1 / median(diff(ts)) / obj.twop.planes.numplanes;

if obj.twop.planes.numplanes > 1
    idx = repmat( (1:obj.twop.planes.numplanes)', [length(ts)/obj.twop.planes.numplanes 1]);
    idx = ismember(idx, obj.twop.planes.planes);
    ts = ts(idx);
    
%     ts = ts(obj.twop.planes.planes: obj.twop.planes.numplanes : end);
%     if length(ts) > size(deconv,1)
%         warning('Number of frames for single plane does not match total number of frames for all planes. Rounding to closest sample size...');
%         ts(end) = [];
%     end
end

obj.twop.ts = ts;

end