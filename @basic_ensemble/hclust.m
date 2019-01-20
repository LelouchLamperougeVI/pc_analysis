function hclust(obj)
% perform agglomerative clustering over SCE and non-SCE epochs to cluster ensembles

% D=sqrt(1-obj.R.^2);
Dm=1-abs(obj.R);
D=squareform(Dm);

obj.tree=linkage(D,'average');

thres=prctile(squareform(1-abs(obj.null_R)),obj.ops.e_prctile);

c=cluster(obj.tree,'cutoff',thres);

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

obj.clust=clust;