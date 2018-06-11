function [I,D,lags] = tdmi(deconv,precision,max_delay)
% Calculates time-delayed mutual information and returns distance matrix

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

edges=[-inf min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
edges(end)=inf;
deconv(isnan(deconv))=min(min(deconv))-1;

%
lags=-max_delay:max_delay;
D=zeros(size(deconv,2),size(deconv,2),2*max_delay+1);
I=zeros(size(deconv,2),size(deconv,2),2*max_delay+1);

parfor i=1:length(lags)
    [I(:,:,i),D(:,:,i)]=get_tdmi(deconv,edges,precision,lags(i));
end


function [I,D]=get_tdmi(deconv,edges,precision,lag)
if lag>0
    deconv1=deconv(1:end-lag,:);
    deconv2=deconv(lag+1:end,:);
else
    deconv1=deconv(-lag+1:end,:);
    deconv2=deconv(1:end+lag,:);
end

p1=zeros(precision+1,size(deconv,2));
p2=zeros(precision+1,size(deconv,2));
p_joint=zeros(size(deconv,2),size(deconv,2),precision+1);
count=1;
for i=1:precision+1
    temp1=deconv1>=edges(i) & deconv1<edges(i+1);
    temp2=deconv2>=edges(i) & deconv2<edges(i+1);
    p1(i,:)=sum(temp1)./size(deconv1,1);
    p2(i,:)=sum(temp2)./size(deconv2,1);
    for j=1:precision+1
        temp2=deconv2>=edges(j) & deconv2<edges(j+1);
        p_joint(:,:,count)=double(temp1')*double(temp2);
        count=count+1;
    end
end
p_joint=p_joint./sum(p_joint,3);
H1=-sum(p1.*log2(p1),'omitnan');
H2=-sum(p2.*log2(p2),'omitnan');


H_joint=-sum(p_joint.*log2(p_joint),3,'omitnan');
I=-H_joint+H1+H2';
I(isnan(I))=0;
D=1-I./H_joint;

I(diag(true(1,length(H1))))=0;
D(diag(true(1,length(H1))))=0;