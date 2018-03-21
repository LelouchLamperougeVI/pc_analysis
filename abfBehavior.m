function behavior=abfBehavior(fn)
% Unfinished code. Temporarily made to calibrate VR and treadmill.

[C,s]=abfload(fn);

frame_ts=C(:,1);
frame_ts=diff(frame_ts)>1;
frame_ts=get_head(frame_ts);
frame_ts=find(frame_ts).*s.*1e-6;

trials=C(:,5);
trials=diff([0;trials])<-.5;
trials=get_head(trials);
trials=find(trials).*s.*1e-6;
trials(trials>frame_ts(end))=[];

trials_ts=arrayfun(@(x) find(frame_ts>=x,1), trials);

cum_pos=behavior.unit_pos;
for i=1:length(behavior.trials)-1
    cum_pos(behavior.trials_ts(i):behavior.trials_ts(i+1)-1)=cum_pos(behavior.trials_ts(i):behavior.trials_ts(i+1)-1)+100*i;
end
cum_pos(behavior.trials_ts(end):end)=cum_pos(behavior.trials_ts(end):end)+100*(i+1);

loc=cum_pos(trials_ts);
gain=mean(diff(loc));

