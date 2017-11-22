function analysis=pc_batch_analysis(behavior,deconv)
% streamlined version of pc_analysis

v2struct(behavior);

vr_length=round(range(unit_pos));
sr=1/mean(diff(frame_ts));

thres=noRun(unit_vel);

sd=4;
bins=50;

[psth,raw_psth,raw_stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,thres,unit_pos,unit_vel,frame_ts,trials);

stack=raw_stack;
stack=(stack-repmat(min(stack),bins,1));
stack=stack./repmat(max(stack),bins,1);

[~,idx]=max(stack);
[~,ordered]=sort(idx);
stack=stack(:,ordered)';

lamb=raw_stack;
m_lamb=mean(lamb);
Pi=Pi./sum(Pi,2);
Pi=Pi';

SI_series=Pi.*lamb./m_lamb.*log2(lamb./m_lamb);
SI=sum(SI_series);

% useGPU=false;
shuffles=1000;
SI=[SI;zeros(shuffles,length(SI))];
parfor i=1:shuffles
    perm=ceil(rand(1)*size(deconv,1));
    shuffled_den=[deconv(perm:end,:);deconv(1:perm-1,:)];
    
    [~,~,lamb1,Pi1]=getStack(bins,sd,vr_length,shuffled_den,thres,behavior.unit_pos,behavior.unit_vel,behavior.frame_ts,behavior.trials);
    m_lamb1=mean(lamb1);
    Pi1=Pi1./sum(Pi1,2);
    Pi1=Pi1';
    
    temp=Pi1.*lamb1./m_lamb1.*log2(lamb1./m_lamb1);
    SI(i+1,:)=sum(temp);
end

pval=1-sum(SI(1,:)>SI(2:end,:))./shuffles;
pc_list=find(pval<0.05);


sparsity=sum(Pi.*lamb).^2./sum(Pi.*lamb.^2);
sparsity=sparsity(1,pc_list);

baseline_thres=range(raw_stack).*.2+min(raw_stack);
width_series=raw_stack>baseline_thres;
width_series=width_series(:,pc_list);

for i=1:size(width_series,2)
    temp=width_series(:,i)';
    start=strfind(temp,[0 1]);
    ending=strfind(temp,[1 0]);
    if temp(1)==1
        start=[1 start];
    end
    if temp(end)==1
        ending=[ending length(temp)];
    end
    temp=ending-start;
    width(i)=max(temp);
end
width=width.*vr_length./bins;

SI=sum(SI_series);
SI=SI(1,pc_list);

analysis=v2struct(vr_length,sr,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width);







