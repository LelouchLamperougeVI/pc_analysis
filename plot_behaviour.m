function h = plot_behaviour(behavior,deconv,analysis)

sig=5;
bins=length(analysis.Pi);

h=figure;
ax1=subplot(4,1,1:2);
deconv=zscore(deconv);
deconv=zscore(fast_smooth(deconv,sig))';
if exist('analysis','var')
    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
    imagesc(-deconv(ordered(end:-1:1),:),'xdata',(behavior.frame_ts-min(behavior.frame_ts))./1);
else
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