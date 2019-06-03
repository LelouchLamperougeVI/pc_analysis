function irCam_ts(obj)
% extract timestamps for IR camera frames

ts=obj.raw(:,obj.get_channel('cam'))>3;
pulse_durations=get_head(ts(end:-1:1));
ts=get_head(ts);
ts=(find(ts) - 1) ./ obj.lfp.fs;
pulse_durations=(find(pulse_durations(end:-1:1)) - 1) ./ obj.lfp.fs;
pulse_durations=pulse_durations-ts;

ts(pulse_durations > 100/obj.lfp.fs)=[];

obj.ts_cam=ts;
