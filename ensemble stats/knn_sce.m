function sce=knn_sce(deconv,i,j)

k=1; %hard limit for now...

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

overlap=false(size(deconv,1),1);
% figure;
% hold on
sce=false(size(deconv,1),1);
pre=-1;I=0;
while(I>pre)
    pre=I;
    
    sce(overlap)=true;
    B=deconv(~sce,:);
    A=B(:,[i j]);
    
    term=sum(sum(isnan(A)).*log(sum(isnan(A))./size(B,1)))-sum(all(isnan(A),2))*log(sum(all(isnan(A),2))/size(B,1));
    
    A=A(all(~isnan(A),2),:); %continuous-continuous space
    
    d=knnsearch(A,A,'k',k+1);
    d=d(:,end); % k th neighbor index
    d=A(d,:); %each row's k-nn
    dist=abs(d-A);
    
    d=abs(B(:,[i j])-permute(A,[3 2 1]));
    n=sum(d<permute(max(dist,[],2),[2 3 1]),1);
    overlap=sum(d<permute(max(dist,[],2),[2 3 1]),3);
    overlap=overlap==max(overlap(~isnan(B(:,[i j])))) & ~isnan(B(:,[i j]));
    overlap=xor(overlap(:,1),overlap(:,2));
    
%     overlap=find(overlap,randi(sum(overlap)));
    
    n=permute(n,[3 2 1]);
    n=n./size(B,1).*size(A,1);
    k_xyz=mean(sum(psi(n+1),2));
    
    I=(psi(k)-k_xyz+psi(size(A,1)))*(size(A,1)/size(B,1))-term/size(B,1);
%     plot(g,I,'o');
end
