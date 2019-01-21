function hclust(obj,varargin)
% perform agglomerative clustering over SCE and non-SCE epochs to cluster ensembles

% D=sqrt(1-obj.R.^2);
obj.set_ops(varargin);
obj.corr;

Dm=1-abs(obj.R);
D=squareform(Dm);

obj.tree=linkage(D,'average');

thres=prctile(squareform(1-abs(obj.null_R)),obj.ops.e_prctile);

c=cluster(obj.tree,'cutoff',thres,'criterion','distance');

count=1;
for i=1:max(c)
    clust{count}=find(c==i)';
    temp=Dm(clust{count},clust{count});
    temp=squareform(temp);
    if mean(temp(:))<thres && length(clust{count})>=obj.ops.e_size
        count=count+1;
    else
        clust(count)=[];
    end
end


% assign each SCE to clusters and remove clusters without SCE
Sm=false(size(obj.deconv,2),length(clust));
for i=1:length(clust)
    Sm(clust{i},i)=true;
end
for i=1:length(obj.SCE.on)
    idx=find(obj.SCE.on(i)==obj.ts):find(obj.SCE.on(i)+obj.SCE.dur(i)==obj.ts);
    mu=obj.deconv(idx,:);
    mu=sum(mu)'.*Sm;
    mu=sum(mu)./cellfun(@length, clust);
    [~,obj.SCE.clust(i)]=max(mu);
end


obj.clust=clust;