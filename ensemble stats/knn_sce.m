function sce=knn_sce(assemblies,deconv,k)

if nargin<3
    k=1;
end

sce=zeros(size(deconv,1),length(assemblies));
for i=1:length(assemblies)
    A=deconv(:,assemblies{i});
    d=knnsearch(A,A,'k',k+1);
    d=d(:,end);
    d=A(d,:); %each row's nn
    dist=abs(d-A);
    d=abs(permute(A,[3 2 1])-A);
    n=sum(d<=dist,3)-1;
    
    loss=sum(psi(n),2);
    [loss,idx]=sort(loss,'descend');
    for j=1:length(loss)
        loss(j)=(size(A,2)-1)*psi(length(loss)-i)-mean(loss(j:end));
    end
    
%     [c,e]=histcounts(loss,'binmethod','fd');
%     [~,thres]=max(c);
%     thres=e(thres+1);
thres=prctile(loss,95);
    loss=loss<thres;
    idx(loss)=[];
    
    sce(idx,i)=1;
end