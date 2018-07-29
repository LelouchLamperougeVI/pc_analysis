function [rho,pval]=fast_spearman(x,y,n,shuffles,gpuFlag)
% A little bit faster than corr(x,y,'type','spearman')
% Expect significant performance increase with GPU

if nargin<4
    shuffles=0;
end
if nargin<5
    gpuFlag=false;
end

if gpuFlag
    A=(A-min(A))./range(A).*(2^16-2);
    A=uint16(A);
    [~,~,ranks]=oddeven(A);
    ranks=double(gather(ranks));
else
    [~,~,x]=par_sort(x);
    [~,~,y]=par_sort(y);
end

rho=zeros(size(x,2),shuffles+1);
rho(:,1)=corr(x,y);
for i=1:shuffles
    perm=repmat(randperm(n)',1,size(x,1)/n);
    perm=perm+(0:size(x,1)/n-1).*n;
    rho(:,i+1)=fast_pearson(x(perm(:),:),y);
end
pval=sum(rho(:,1)>rho(:,2:end),2)./shuffles;
rho=rho(:,1);