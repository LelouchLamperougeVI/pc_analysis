function two_photon_ts(obj, deconv)
% extract 2p timestamps

ts = obj.abf.raw(:,obj.get_channel('2p')) < 3;
ts = get_head(ts);
ts = (find(ts) - 1)./obj.lfp.fs;
if ts(1) == 0
    ts(1) = [];
end

% For future me, this here below is a better approach for determining frame
% times. The scanner spends ~50 ms to scan through a frame. Right now we
% are counting the beginning of the scan as the frame timestamps. But if we
% take the middle of the scan, we get true $\pm$ 25ms confidence interval
% for ~20HZ sampling rate. So far we haven't encountered a problem
% requiring comparing LFP and imaging times directly. But in the future
% please remember to turn this on and make the same changes where
% appropriate (easy to identify using grep -ri "get_channel('2p')").

% ts = obj.abf.raw(:,obj.get_channel('2p')) > 4;
% heads = get_head(ts);
% tails = get_head(ts(end:-1:1)); tails = tails(end:-1:1);
% ts = (find(heads) + ((find(tails) - find(heads)) ./ 2)) ./ obj.lfp.fs;

obj.twop.fs = 1 / median(diff(ts)) / obj.twop.planes.numplanes;

if obj.twop.planes.numplanes > 1
%     idx = repmat( (1:obj.twop.planes.numplanes)', [length(ts)/obj.twop.planes.numplanes 1]);
    idx = repmat( (1:obj.twop.planes.numplanes)', [size(deconv, 1)/obj.twop.planes.numplanes 1]);
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