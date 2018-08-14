function I=knn_mi(deconv,k)
% Implementation of k-nearest neighbor based MI estimator

if nargin<2
    k=1;
end

[~,~,deconv]=par_sort(deconv); %rank-order data
% deconv=ca_filt(deconv);
% deconv=log(deconv);
% deconv(isinf(deconv))=nan;
% deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan')+100;
% deconv(isnan(deconv))=0;

h=waitbar(0,'estimating MI using k-NN estimator...');
I=zeros(size(deconv,2));
for i=1:size(deconv,2)
    for j=i+1:size(deconv,2)
        A=deconv(:,[i j]);
        A(sum(A,2)==0,:)=[];
        d=knnsearch(A,A,'k',k+1);
        d=d(:,end);
        d=A(d,:); %each row's k-nn
%         dist=abs(d-A);
        dist=max(abs(d-A),[],2);
        d=abs(permute(A,[3 2 1])-A);
        n=sum(d<dist,3);
        
%         I(i,j)=psi(k)-1/k-mean(sum(psi(n),2))+psi(size(deconv,1));
        I(i,j)=psi(k)-mean(sum(psi(n+1),2))+psi(size(A,1));
    end
    waitbar(i/size(deconv,2),h);
end
close(h);

for i=1:size(deconv,2)
    for j=i+1:size(deconv,2)
        I(j,i)=I(i,j);
    end
end