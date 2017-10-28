function [psth,raw_psth,raw_stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,vel_thres,unit_pos,unit_vel,frame_ts,trials)

sd=sd/vr_length*bins;
list=1:size(deconv,2);

vel_thres=find(unit_vel>vel_thres | unit_vel<-vel_thres);
vel=unit_vel(vel_thres);
pos=unit_pos(vel_thres);
ft=frame_ts(vel_thres);

% raw_psth=arrayfun(@(x) zeros(length(trials)-1,bins),1:length(list),'uniformoutput',false);
raw_psth=zeros(length(trials)-1,bins,length(list));
vel_stack=zeros(length(trials)-1,bins);

edges=zeros(1,(bins+1)*(length(trials)-1));
idx=min(pos):range(pos)/(bins+1):max(pos);
count=1;
for i=1:length(trials)-1
    for j=1:bins+1
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
end
occupancy_series=diff(edges);
occupancy_series=repmat(occupancy_series,length(list),1);
Pi=zeros(size(deconv,2),bins); %for SI test

% for k=1:length(list)
%     n=list(k);
    signal=deconv(vel_thres,:);
        
    count1=1;
    for i=1:length(trials)-1
        for j=1:bins
            raw_psth(i,j,:)=reshape(mean(signal(edges(count1):edges(count1+1),:),1),1,1,length(list));
            raw_psth(i,j,:)=raw_psth(i,j,:)./reshape(occupancy_series(:,count1),1,1,length(list));
            idx=reshape(isinf(raw_psth(i,j,:)),1,length(list));
            
            occupancy_series(idx,count1)=1;
            raw_psth(i,j,:)=reshape(mean(signal(edges(count1):edges(count1+1),:),1),1,1,length(list));
            raw_psth(i,j,:)=raw_psth(i,j,:)./reshape(occupancy_series(:,count1),1,1,length(list));
            
            vel_stack(i,j)=mean(vel(edges(count1):edges(count1+1)));
            
            Pi(:,j)=occupancy_series(:,count1);
            count1=count1+1;
        end
    count1=count1+1;
    end
%     raw_psth{k}(all(isnan(raw_psth{k}),2),:)=[];
    raw_psth(isnan(raw_psth))=0;
    psth=arrayfun(@(x) Smooth(raw_psth(:,:,x),[0 sd]),1:length(list),'uniformoutput',false);
% end
% Pi=occupancy_series;
stack=arrayfun(@(x) mean(raw_psth(:,:,x)),1:length(list),'uniformoutput',false);
stack=cell2mat(stack);
stack=reshape(stack,bins,size(deconv,2));
stack=Smooth(stack,[sd 0]);
raw_stack=stack;