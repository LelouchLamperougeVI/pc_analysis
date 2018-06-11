function [I,D]=get_mi(deconv,precision)
% returns MI and distance dissimilarity matrices
% D = 1 - I(x;y) / H(x,y)
%
% Now works with both binary and Q matrices
% p = {P(x=0), P(precision<log(x)<precision), ...}
% 
% if ~islogical(deconv)
%     deconv=logical(deconv);
% end
if islogical(deconv)
    p=sum(deconv)./size(deconv,1);
    H=-(p.*log2(p)+(1-p).*log2(1-p));

    p_joint=zeros(size(deconv,2),size(deconv,2),4);

    p_joint(:,:,1)=double(deconv')*double(deconv);
    p_joint(:,:,2)=double(~deconv')*double(deconv);
    p_joint(:,:,3)=double(deconv')*double(~deconv);
    p_joint(:,:,4)=double(~deconv')*double(~deconv);
    p_joint=p_joint./sum(p_joint,3);
else
    if ~exist('precision','var')
        error('precision required for Q-matrix type input');
    end
    
    deconv=ca_filt(deconv);
    deconv=log(deconv);
    deconv(isinf(deconv))=nan;
    deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');
    
    edges=[-inf min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
    edges(end)=inf;
    deconv(isnan(deconv))=min(min(deconv))-1;
    
    p=zeros(precision+1,size(deconv,2));
    p_joint=zeros(size(deconv,2),size(deconv,2),(precision+1)^2);
    count=1;
    for i=1:precision+1
        temp1=deconv>=edges(i) & deconv<edges(i+1);
        p(i,:)=sum(temp1)./size(deconv,1);
        for j=1:precision+1
            temp2=deconv>=edges(j) & deconv<edges(j+1);
            p_joint(:,:,count)=double(temp1')*double(temp2);
            count=count+1;
        end
    end
    p_joint=p_joint./sum(p_joint,3);
    H=-sum(p.*log2(p),'omitnan');
end

H_joint=-sum(p_joint.*log2(p_joint),3,'omitnan');
I=-H_joint+H+H';
I(isnan(I))=0;
D=1-I./H_joint;

I(diag(true(1,length(H))))=0;
D(diag(true(1,length(H))))=0;