function [psth,raw_psth,raw_count,raw_stack,stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,vel_thres,unit_pos,unit_vel,frame_ts,trials)

sd=sd/vr_length*bins;
list=1:size(deconv,2);

vel_thres=unit_vel>vel_thres | unit_vel<-vel_thres;
vel=unit_vel(vel_thres);
pos=unit_pos(vel_thres);
ft=frame_ts(vel_thres);

% raw_psth=arrayfun(@(x) zeros(length(trials)-1,bins),1:length(list),'uniformoutput',false);
raw_psth=zeros(length(trials)-1,bins,length(list));
raw_count=raw_psth;
vel_stack=zeros(length(trials)-1,bins);

edges=zeros(1,(bins+1)*(length(trials)-1));
idx=min(pos):range(pos)/(bins):max(pos);
count=1;
for i=1:length(trials)-1
    for j=1:bins
        try
            edges(count)=find(pos>idx(j) & pos<=idx(j+1) & ft >= trials(i) & ft <= trials(i+1),1);
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
signal=ca_filt(signal);
signal=double(signal>0);
count1=1;
for i=1:length(trials)-1
    for j=1:bins
%         temp=mean(signal(edges(count1):edges(count1+1),:),1);
        temp=sum(signal(edges(count1):edges(count1+1),:),1);
        raw_count(i,j,:)=temp;
        raw_psth(i,j,:)=temp./occupancy_series(count1);
%         if any(isinf(raw_psth(i,j,:)))
%             error(['fatal error: 0 occupancy in trial ' num2str(i)]);
%         end

        vel_stack(i,j)=mean(vel(edges(count1):edges(count1+1)));

        Pi(j)=Pi(j)+occupancy_series(count1);
        count1=count1+1;
    end
    count1=count1+1;
end
raw_psth(isnan(raw_psth) | isinf(raw_psth))=0;
% Pi=Pi./sum(Pi);

psth=arrayfun(@(x) fast_smooth(raw_psth(:,:,x),sd,2),1:length(list),'uniformoutput',false);

raw_stack=arrayfun(@(x) mean(raw_psth(:,:,x)),1:length(list),'uniformoutput',false);
raw_stack=cell2mat(raw_stack);
raw_stack=reshape(raw_stack,bins,size(deconv,2));

stack=arrayfun(@(x) mean(fast_smooth(raw_psth(:,:,x),sd,2)),1:length(list),'uniformoutput',false);
stack=cell2mat(stack);
stack=reshape(stack,bins,size(deconv,2));
stack=(stack-min(stack));
stack=stack./max(stack);