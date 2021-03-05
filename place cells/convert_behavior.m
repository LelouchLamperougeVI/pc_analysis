function [behavior,deconv]=convert_behavior(behavior,tcs,deconv)
%convert Dun's analysis into mine

frame_ts=tcs.tt;
idx=find(frame_ts>behavior.ts(end),1);
frame_ts(idx:end)=[];
deconv(idx:end,:)=[];
idx=find(frame_ts<behavior.ts(1),1,'last');
frame_ts(1:idx)=[];
deconv(1:idx,:)=[];

trials_ts=arrayfun(@(x) find(behavior.trial==x,1),1:max(behavior.trial));
trials_ts=behavior.ts(trials_ts);
trials_ts(trials_ts>frame_ts(end))=[];
trials_ts(trials_ts<frame_ts(1))=[];
% trials_ts=arrayfun(@(x) find(frame_ts>=x,1),trials_ts);
trials_ts=knnsearch(frame_ts',trials_ts');
trials=frame_ts(trials_ts);

% idx=arrayfun(@(x) find(behavior.ts>=x,1),frame_ts);
idx=knnsearch(behavior.ts',frame_ts');

unit_pos=behavior.pos_raw(idx);
% unit_vel=behavior.speed(idx);
unit_vel = diff([0 behavior.pos_cum(idx)])./ median(diff(frame_ts));

clear behavior;

behavior.frame_ts=frame_ts;
behavior.trials=trials;
behavior.trials_ts=trials_ts;
behavior.unit_pos=unit_pos;
behavior.unit_vel=unit_vel;