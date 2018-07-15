function [psth,raw_psth,raw_count,edges,raw_stack,stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,vel_thres,unit_pos,unit_vel,frame_ts,trials)
% looking back at my old code from a year ago, it's such a sloppy
% inneficient POS...

sd=sd/vr_length*bins;
list=1:size(deconv,2);

vel_thres=unit_vel>vel_thres | unit_vel<-vel_thres;
vel=unit_vel(vel_thres);
pos=unit_pos(vel_thres);
ft=frame_ts(vel_thres);

raw_psth=zeros(length(trials)-1,bins,length(list));
raw_count=cell(1,bins);
vel_stack=zeros(length(trials)-1,bins);

edges=zeros(1,(bins+1)*(length(trials)-1));
idx=min(pos):range(pos)/(bins):max(pos);
count=1;
for i=1:length(trials)-1
    for j=1:bins
        try
            temp=find(pos>idx(j) & pos<=idx(j+1) & ft >= trials(i) & ft <= trials(i+1),1);
            edges(count)=temp(get_head(temp>edges(count-1)));
        catch
            try
                edges(count)=edges(count-1);
            catch
                edges(count)=find(ft>=trials(1),1);
            end
        end
        count=count+1;
    end
    [~,edges(count)]=min(abs(ft-trials(i+1)));
    count=count+1;
end
occupancy_series=diff(edges);
Pi=zeros(1,bins); %for SI test

signal=deconv(vel_thres,:);
signal_log=ca_filt(signal);
signal_log=log(signal_log);
signal_log(isinf(signal_log))=nan;
signal_log=(signal_log-mean(signal_log,'omitnan'))./std(signal_log,'omitnan');

raw_stack=zeros(bins,size(deconv,2));
% raw_stack_norm=raw_stack;
for i=1:length(idx)-1
    temp=pos>idx(i) & pos<=idx(i+1) & ft >= trials(1) & ft <= trials(end);
    raw_stack(i,:)=mean(signal(temp,:));
%     raw_stack_norm(i,:)=mean(signal(temp,:))./sum(temp); %normalize by occupancy
end
raw_stack(isnan(raw_stack))=0;
% raw_stack_norm(isnan(raw_stack_norm))=0;

count1=1;
for i=1:length(trials)-1
    for j=1:bins
%         temp=mean(signal(edges(count1):edges(count1+1),:),1);
%         temp=signal(edges(count1):edges(count1+1),:);
%         raw_count(i,j,:)=temp;
        raw_count{j}=[raw_count{j}; signal_log(edges(count1):edges(count1+1)-1,:)];
%         raw_psth(i,j,:)=mean(signal(edges(count1):edges(count1+1),:),1)./(occupancy_series(count1).^(1/3)); % cube root transformation to reduce strong positive skew in occupancy distribution
        raw_psth(i,j,:)=mean(signal(edges(count1):edges(count1+1)-1,:),1);

        vel_stack(i,j)=mean(vel(edges(count1):edges(count1+1)-1));

        Pi(j)=Pi(j)+occupancy_series(count1);
        count1=count1+1;
    end
    count1=count1+1;
end
nbins=floor(mean(fd_bins(signal_log))); %using Freedman-Diaconis's rule to determine optimal # of bins
edges=linspace(min(signal_log(:)),max(signal_log(:)),nbins);
raw_psth(isnan(raw_psth) | isinf(raw_psth))=0;

psth=arrayfun(@(x) fast_smooth(raw_psth(:,:,x),sd,2),list,'uniformoutput',false);

% raw_stack=mean(raw_psth)./Pi;
% raw_stack=squeeze(raw_stack);
% raw_stack=raw_stack./Pi';
% raw_psth=raw_psth_norm;

% stack=arrayfun(@(x) mean(fast_smooth(raw_psth(:,:,x),sd,2)),list,'uniformoutput',false);
% stack=cell2mat(stack);
% stack=reshape(stack,bins,size(deconv,2));
stack=fast_smooth(raw_stack,sd);
% stack=fast_smooth(raw_stack_norm,sd);
% stack=raw_stack;
stack=(stack-min(stack));
stack=stack./max(stack);