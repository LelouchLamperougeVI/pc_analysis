c=3;

%% run
dec=fast_smooth(deconv(:,ass.clust{c}),6);
r=xcorr(dec,30);
r=r(:,triang_idx);
% [~,x]=max(r); [~,idx]=sort(x);
x=max(r); [~,idx]=sort(x);
figure
% imagesc(fast_smooth(zscore(r(:,idx)),0)')
imagesc(fast_smooth(r(:,idx),0)')

[~,x]=max(fast_smooth(zscore(r(:,idx)),0));
hold on
plot(x,1:length(idx),'r*');

%% rest
dec=fast_smooth(deconv(:,ass.clust{c}),6);
r=xcorr(dec,30);
r=r(:,triang_idx);
figure
% imagesc(fast_smooth(zscore(r(:,idx)),0)')
imagesc(fast_smooth(r(:,idx),0)')

[~,x]=max(fast_smooth(zscore(r(:,idx)),0));
hold on
plot(x,1:length(idx),'r*');

%% up triangle idx
N=size(dec,2);
N=reshape(1:N^2,N,N);
triang_idx=tril(N,-1);
triang_idx=triang_idx(triang_idx>0);

%%
order=get_order(analysis);
order=intersect(order,ass.clust{c},'stable');
figure
imagesc(fast_smooth(zscore(deconv(:,order)),5)');
colormap jet


%% average sequence
temp=0;
for i = 1:length(ass.clust_SCE(c).SCE.peak)
    idx = find(ass.clust_SCE(c).SCE.peak(i) == ass.ts);
    temp = temp + deconv(idx-30:idx+30,:);
end

%% sequences side-by-side
temp=[];
order=get_order(analysis);
for i = 1:length(ass.clust_SCE(c).SCE.peak)
    idx = find(ass.clust_SCE(c).SCE.peak(i) == ass.ts);
    
%     wdw = deconv(idx-30:idx+30,intersect(order,ass.clust{c},'stable'));
%     [~,idx]=max(fast_smooth(wdw,2));
%     idx = sub2ind(size(wdw),idx,1:size(wdw,2));
%     wdw=zeros(size(wdw));
%     wdw(idx)=1;
%     temp = [temp; ones(1,size(wdw,2)); wdw];
    
    temp = [temp; deconv(idx-30:idx+30,intersect(order,ass.clust{c},'stable'))];
end

