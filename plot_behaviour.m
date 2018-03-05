function h = plot_behaviour(behavior,deconv,analysis)

sig=5;

h=figure;
ax1=subplot(4,1,1:2);
deconv=zscore(fast_smooth(deconv,sig))';
if exist('analysis','var')
    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
    imagesc(deconv(:,ordered));
else
    imagesc(deconv);
end
ylabel('neuron no.');
ax2=subplot(4,1,3);
plot(behavior.unit_pos);
ylabel('position (cm)');
ax3=subplot(4,1,4);
plot(behavior.unit_vel);
xlabel('frames');
ylabel('velocity (cm/s');

linkaxes([ax1,ax2,ax3],'x');