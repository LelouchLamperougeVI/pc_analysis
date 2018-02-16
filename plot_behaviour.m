function h = plot_behaviour(behavior,deconv,sig)
if nargin==2
    sig=5;
end
h=figure;
ax1=subplot(4,1,1:2);
imagesc(fast_smooth(deconv,sig)');
ylabel('neuron no.');
ax2=subplot(4,1,3);
plot(behavior.unit_pos);
ylabel('position (cm)');
ax3=subplot(4,1,4);
plot(behavior.unit_vel);
xlabel('frames');
ylabel('velocity (cm/s');

linkaxes([ax1,ax2,ax3],'x');