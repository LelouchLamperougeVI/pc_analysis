function two_photon_ts(obj)
% extract 2p timestamps

ts = obj.raw(:,obj.get_channel('2p')) < 3;
ts = get_head(ts);
ts = (find(ts) - 1)./obj.lfp.fs;
if ts(1) == 0
    ts(1) = [];
end
obj.twop.ts = ts;
obj.twop.fs = 1/median(diff(ts));

end