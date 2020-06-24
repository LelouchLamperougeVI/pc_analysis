function hclust(obj,varargin)
% perform agglomerative clustering over SCE and non-SCE epochs to cluster ensembles

% D=sqrt(1-obj.ensembles.R.^2);
obj.set_ops(varargin);
if isempty(obj.ensembles.R)
    obj.corr;
end

% Dm=1-abs(obj.ensembles.R); %this used to be default
Dm = 1 - obj.ensembles.R; %this is the default when using 'correlation' for pdist
Dm(1:length(Dm)+1:numel(Dm))=0;
D=squareform(Dm);

obj.ensembles.tree=linkage(D,'average');
try
    obj.ensembles.clust_order = optimalleaforder(obj.ensembles.tree, D);
catch
    tmp_d=D; tmp_d(isnan(tmp_d)) = 1;
    obj.ensembles.clust_order = optimalleaforder(obj.ensembles.tree, tmp_d);
end

try
    thres = linkage(squareform(1-abs(obj.ensembles.null_R)), 'average');
    thres=prctile(thres(:,3),obj.ops.e_prctile);
    obj.ensembles.h_thres = thres;
catch
end

if strcmp(obj.ops.clust_method, 'silhouette')
    c=obj.silhouette_cluster(obj.ensembles.tree, Dm, obj.ops.e_size);
elseif strcmp(obj.ops.clust_method, 'shuffle')
    c=cluster(obj.ensembles.tree,'cutoff',thres,'criterion','distance');
elseif strcmp(obj.ops.clust_method, 'thres')
    c=cluster(obj.ensembles.tree,'cutoff',obj.ops.clust_thres,'criterion','distance');
end

count=1;
for i=1:max(c)
    clust{count}=find(c==i)';
    temp=Dm(clust{count},clust{count});
    temp=squareform(temp);
%     if mean(temp(:))<thres && length(clust{count})>=obj.ops.e_size
    if length(clust{count})>=obj.ops.e_size
        count=count+1;
    else
        clust(count)=[];
    end
end


% assign each SCE to clusters and remove clusters without SCE
if isempty(clust)
    obj.ensembles.clust=cell(0);
    return
end

Sm=false(size(obj.twop.deconv,2),length(clust));
for i=1:length(clust)
    Sm(clust{i},i)=true;
end
if ~isempty(obj.ensembles.SCE)
    for i=1:length(obj.ensembles.SCE.on)
        idx=find(obj.ensembles.SCE.on(i)==obj.twop.ts):find(obj.ensembles.SCE.on(i)+obj.ensembles.SCE.dur(i)==obj.twop.ts);
        mu=obj.twop.deconv(idx,:);
        mu=sum(mu)'.*Sm;
        mu=sum(mu)./cellfun(@length, clust);
        [~,obj.ensembles.SCE.clust(i)]=max(mu);
    end
end

obj.ensembles.clust=clust;
obj.make_colours;

if strcmpi(obj.ops.order,'cluster')
    obj.set_ops('order','cluster');
end