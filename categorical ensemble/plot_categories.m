function h=plot_categories(behavior,deconv,sig)
if nargin<3
    sig=10;
end
h=figure;
ax1=subplot(2,1,1);
imagesc(fast_smooth(deconv,sig)');
ax2=subplot(2,1,2);
plot(behavior.object_frame);
hold on
plot(behavior.black_frame);
legend('objects onset','blank onset');
linkaxes([ax1 ax2],'x');

idx=find(behavior.object_frame);
for i=1:length(behavior.object_ts)
    text(idx(i),1-i/length(idx),behavior.object_type(i,:));
end
