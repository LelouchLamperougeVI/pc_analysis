function h = plot_behaviour(analysis,deconv,peaks_flag)

if nargin < 3
    peaks_flag = false;
end

sig=3;

h=figure;
ax1=subplot(4,1,1:2);

bins=length(analysis.Pi);
behavior=analysis.behavior;

run_thres = noRun(behavior.unit_vel);
run_thres = behavior.unit_vel < run_thres;

deconv=zscore(fast_smooth(deconv,sig));

ordered = get_order(analysis);
ordered = intersect(ordered, analysis.pc_list, 'stable');

if peaks_flag
    p_neurons = repmat(1:length(ordered), length(behavior.trials)-1);
    p_frame = zeros(1, length(ordered) * (length(behavior.trials) - 1));
    for i = 1:length(behavior.trials)-1
        temp = deconv(:, ordered);
        temp(run_thres, :) = nan;
        temp(setxor(1:size(deconv,1), behavior.trials_ts(i):behavior.trials_ts(i+1)-1), :) = nan;
        [~,idx] = max(temp);
        p_frame((i-1)*length(ordered)+1 : i*length(ordered)) = idx;
    end
    p_frame = behavior.frame_ts(p_frame)-min(behavior.frame_ts);
    plot(p_frame, p_neurons, 'k.', 'markersize',20);
    ylabel('neuron no.');
else
    imagesc(-deconv(:, ordered(end:-1:1))','xdata',(behavior.frame_ts-min(behavior.frame_ts))./1);
    colormap gray
    ylabel('neuron no.');
end

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
% if ~exist('range','var')
idx=linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins+1);
% else
%     idx=linspace(-100,0,bins+1);
% end
for i=1:length(behavior.trials)-1
    for j=1:length(idx)-1
        tmp=behavior.unit_pos>idx(j) & behavior.unit_pos<=idx(j+1) & behavior.frame_ts>behavior.trials(i) & behavior.frame_ts<=behavior.trials(i+1);
        vel(j,i)=mean(behavior.unit_vel(tmp));
    end
end
sd=4/analysis.vr_length*bins;
for i=1:size(vel,2)
    if isnan(vel(end,i))
        idx=get_head(isnan(vel(:,i)));
        idx=find(idx,1,'last');
        vel(idx:end,i)=0;
    end
    if isnan(vel(1,i))
        idx=get_head(flipud(isnan(vel(:,i))));
        idx=flipud(idx);
        idx=find(idx,1);
        vel(1:idx,i)=0;
    end
end
vel=fillmissing(vel,'linear');
vel=fast_smooth(vel,sd);
plot(linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins),vel,'color',[.5 .5 .5]);
hold on
plot(linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins),mean(vel,2),'k','linewidth',2);
xlabel('position (cm)');
ylabel('velocity (cm/s)');
xlim([-100 0])
ylim([0 50])
