function [sce,loss]=mi_sce(assemblies,deconv,max_it,jitter)
if nargin<3
    max_it=10;
end
if nargin<4
    jitter=10;
end

sce=zeros(size(deconv,1),length(assemblies));
% idx=zeros(size(deconv,1),length(assemblies));
loss=cell(1,length(assemblies));
for i=1:length(assemblies)
    temp=deconv(:,assemblies{i});
    [~,cutoff]=get_mi(temp);
    cutoff=triu(cutoff,1);
    cutoff=mean(cutoff(cutoff~=0));
    
    [~,idx]=sort(sum(temp,2),'descend');
    temp=temp(idx,:);
    
    cost=0;
    it=0;
    count=1;
    while it<=max_it && count<=size(temp,1)
%     while cost(end)<cutoff && it<=max_it && count<=size(temp,1)
        [~,d]=get_mi(temp(1:count,:));
        d=triu(d,1);
        d=mean(d(d~=0));
        d(isnan(d))=-inf;
        cost=[cost d];
        if cost(end)<=cost(end-1)
            it=it+1;
        end
        count=count+1;
    end
    cost(1)=[];
    loss{i}=cost;
    
    sce(idx(1:count-1),i)=1;
    idx=find(sce(:,i));
    idx2=find(diff(idx)<=jitter);
    for j=1:length(idx2)
        sce(idx(idx2(j)):idx(idx2(j)+1),i)=1;
    end
    
    figure;
    ax1=subplot(2,1,1);
    imagesc(deconv(:,assemblies{i})');
    ax2=subplot(2,1,2);
    plot(sce(:,i));
    linkaxes([ax1 ax2],'x');
end