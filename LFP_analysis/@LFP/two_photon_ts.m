function two_photon_ts(obj, deconv)
% extract 2p timestamps

ts = obj.raw(:,obj.get_channel('2p')) < 3;
ts = get_head(ts);
ts = (find(ts) - 1)./obj.lfp.fs;
if ts(1) == 0
    ts(1) = [];
end

if obj.twop.numplanes > 1
    ts = ts(obj.twop.plane: obj.twop.numplanes : end);
    if length(ts) > size(deconv,1)
        warning('Number of frames for single plane does not match total number of frames for all planes. Rounding to closest sample size...');
        ts(end) = [];
    end
end

obj.twop.ts = ts;
obj.twop.fs = 1/median(diff(ts));

end