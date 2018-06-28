function h = plot_behaviour(behavior,deconv)

sig=5;
bins=80;

h=figure;
ax1=subplot(4,1,1:2);
if nargin<2
    analysis=behavior;
    bins=length(analysis.Pi);
    behavior=analysis.behavior;
    deconv=analysis.deconv;
    
    deconv=zscore(deconv);
    deconv=zscore(fast_smooth(deconv,sig))';

    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
    imagesc(-deconv(ordered(end:-1:1),:),'xdata',(behavior.frame_ts-min(behavior.frame_ts))./1);
else
    deconv=zscore(deconv);
    deconv=zscore(fast_smooth(deconv,sig))';
    
    imagesc(-deconv,'xdata',(behavior.frame_ts-min(behavior.frame_ts))./1);
end
colormap gray
ylabel('neuron no.');
ax2=subplot(4,1,3);
plot((behavior.frame_ts-min(behavior.frame_ts))./1,behavior.unit_pos);
ylabel('position (cm)');
ax3=subplot(4,1,4);
plot((behavior.frame_ts-min(behavior.frame_ts))./1,behavior.unit_vel);
xlabel('time (sec)');
ylabel('velocity (cm/s');

linkaxes([ax1,ax2,ax3],'x');

figure;
vel_thres=noRun(behavior.unit_vel);
vel_thres=behavior.unit_vel>vel_thres | behavior.unit_vel<-vel_thres;
behavior.unit_vel=behavior.unit_vel(vel_thres);
behavior.unit_pos=behavior.unit_pos(vel_thres);
behavior.frame_ts=behavior.frame_ts(vel_thres);

vel=zeros(bins,length(behavior.trials)-1);
idx=linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins+1);
for i=1:length(behavior.trials)-1
    for j=1:length(idx)-1
        tmp=behavior.unit_pos>idx(j) & behavior.unit_pos<=idx(j+1) & behavior.frame_ts>behavior.trials(i) & behavior.frame_ts<=behavior.trials(i+1);
        vel(j,i)=mean(behavior.unit_vel(tmp));
    end
end
sd=4/analysis.vr_length*bins;
vel=fillmissing(vel,'linear');
vel=fast_smooth(vel,sd);
plot(linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins),vel,'color',[.5 .5 .5]);
hold on
plot(linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins),mean(vel,2),'k','linewidth',2);
xlabel('position (cm)');
ylabel('velocity (cm/s)');
xlim([-100 0])
ylim([0 50])
