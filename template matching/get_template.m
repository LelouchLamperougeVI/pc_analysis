function template=get_template(analysis)

sd=0;

behavior=analysis.behavior;
deconv=analysis.deconv;

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

vr_length=round(range(unit_pos));
thres=noRun(unit_vel);

thres=unit_vel>thres | unit_vel<-thres;
ft=frame_ts(thres);
bins=knnsearch(ft',trials');

bins=round(median(diff(bins)));

[~,~,~,~,template]=getStack(bins,sd,vr_length,deconv,thres,unit_pos,unit_vel,frame_ts,trials);