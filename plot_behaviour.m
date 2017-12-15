function h = plot_behaviour(behavior,deconv,sig)
if nargin==2
    sig=5;
end
h=figure;
subplot(4,1,1:2);
imagesc(fast_smooth(deconv,sig)');
ylabel('neuron no.');
subplot(4,1,3);
plot(behavior.unit_pos);
ylabel('position (cm)');
subplot(4,1,4);
plot(behavior.unit_vel);
xlabel('frames');
ylabel('velocity (cm/s');