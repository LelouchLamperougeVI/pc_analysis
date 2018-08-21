function [I,D]=knn_mi(deconv,k,thres)
% Implementation of k-nearest neighbor based MI estimator specifically for
% calcium data

if nargin<2
    k=1;
end
if nargin<3
    thres=k+1;
end

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

lol=figure;

h=waitbar(0,'estimating MI using k-NN estimator...');
I=NaN(size(deconv,2));
D=I;
for i=1:size(deconv,2)
    for j=i+1:size(deconv,2)
        A=deconv(:,[i j]);
        
        term=sum(sum(isnan(A)).*log(sum(isnan(A))./size(deconv,1)))-sum(all(isnan(A),2))*log(sum(all(isnan(A),2))/size(deconv,1));
        n_xy=sum(~isnan(A));
        n_xy0=sum(isnan(A));
        
        A=A(all(~isnan(A),2),:);
        if size(A,1)>=thres
            d=knnsearch(A,A,'k',k+1);
            d=d(:,end);
            d=A(d,:); %each row's k-nn
            dist=abs(d-A);
%             d=abs(permute(A,[3 2 1])-A);
            d=abs(deconv(:,[i j])-permute(A,[3 2 1]));
            n=sum(d<permute(max(dist,[],2),[2 3 1]),1);
            n=permute(n,[3 2 1]);
            
            k_xy=mean(psi(n+1));
            
            n=(n+1)./size(deconv,1).*size(A,1);
            k_xyz=mean(sum(psi(n+1),2));
            
            H_xy=n_xy./size(deconv,1).*(-k_xy+psi(n_xy)+mean(dist)./n_xy)-n_xy0./size(deconv,1).*log(n_xy0./size(deconv,1));
            
            I(i,j)=(psi(k)-k_xyz+psi(size(A,1)))*(size(A,1)/size(deconv,1))-term/size(deconv,1);
%             I(i,j)=(psi(k)-k_xy+psi(size(A,1)))*(size(A,1)/size(deconv,1));
            D(i,j)=1-I(i,j)/max(H_xy);
        else
            I(i,j)=nan;
        end
    end
    figure(lol);
    imagesc(I);
    waitbar(i/size(deconv,2),h);
end
close(h);
close(lol);

for i=1:size(deconv,2)
    for j=i+1:size(deconv,2)
        I(j,i)=I(i,j);
        D(j,i)=D(i,j);
    end
end

D(isnan(D))=1;
for i=1:size(deconv,2)
    D(i,i)=0;
end