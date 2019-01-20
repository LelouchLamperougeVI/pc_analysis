function plot(obj)
%Plot results from detection

figure;
if ~isempty(obj.order)
    order=obj.order;
    lab='sorted neuron no.';
else
    order=1:size(obj.deconv,2);
    lab='neuron no.';
end
deconv=fast_smooth(obj.deconv(:,order),obj.ops.sig);
ax(1)=subplot(4,1,1:3);
imagesc('xdata',obj.ts,'cdata',deconv');
ylim([1 size(deconv,2)]);
ylabel(lab);

ax(2)=subplot(4,1,4);
plot(obj.ts,obj.MUA);
hold on
plot(obj.SCE.on, min(obj.MUA).*ones(length(obj.SCE.on),1),'^k');
ylabel('population mean dF/F');
xlabel('time (sec)');

linkaxes(ax,'x');