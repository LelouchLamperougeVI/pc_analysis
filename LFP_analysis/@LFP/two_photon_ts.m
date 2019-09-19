function two_photon_ts(obj, deconv)
% extract 2p timestamps

ts = obj.raw(:,obj.get_channel('2p')) < 3;
ts = get_head(ts);
ts = (find(ts) - 1)./obj.lfp.fs;
if ts(1) == 0
    ts(1) = [];
end

if obj.twop.numplanes > 1
    if mod(length(ts), size(deconv,1))
        error('number of frames for single plane does not match total number of frames for all planes');
    end
    ts = ts(obj.twop.plane: obj.twop.numplanes : end);
end

obj.twop.ts = ts;
obj.twop.fs = 1/median(diff(ts));

end