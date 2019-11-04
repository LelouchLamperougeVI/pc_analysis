delay = 500;

tcs.tt=lfp.twop.ts';
[behavior,deconv]=convert_behavior(lfp.behavior,tcs,lfp.twop.deconv);
bins = 50;
sd = 4;

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

thres=noRun(unit_vel);
thres=(unit_vel>thres | unit_vel<-thres) & (trials(1) < frame_ts & trials(end) > frame_ts);
unit_vel=unit_vel(thres);
unit_pos=unit_pos(thres);
frame_ts=frame_ts(thres);

deconv=ca_filt(deconv);
deconv=deconv(thres,:);
vr_length=round(range(unit_pos));

SI = zeros(size(deconv,2), delay*2+1);
tic
parfor ii = 1:delay*2+1
    [~,~,raw_stack,mu_fr,Pi]=getStack(bins,sd,vr_length, circshift(deconv, ii-delay-1, 1) ,unit_pos,unit_vel,frame_ts,trials);
    SI(:, ii)=get_si_skaggs(raw_stack,mu_fr,Pi);
end
toc

%%
ratio = sum( SI(:,delay+1) > SI) ./ size(SI,1);
% ratio = sum( SI(lfp.analysis.pc_list,delay+1) > SI(lfp.analysis.pc_list,:)) ./ length(lfp.analysis.pc_list);
% t = (-delay:delay) ./ lfp.twop.fs;
t = (-delay:delay);
figure; hold on
plot(t, ratio, 'k');
xlim([-600 600]);
ylim([0 1]);
area(-delay:-251, ratio(1:250), 'facecolor','blue');
area(251:delay, ratio(end-249:end), 'facecolor','blue');
area(-250:250, ratio(250:750), 'facecolor','red');

xlabel('time delay (frames)');
ylabel('fraction SI_0 > SI_t');

%%
Aa = sum(ratio([1:500 502:end]));
Ab = sum(ratio([1:250 end-250+1:end]));
bias = [bias; 1 - ( Aa + Ab*11000/500 ) / ( Ab / 500 * 12000)];