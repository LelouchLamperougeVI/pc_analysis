function [I,D]=get_mi(deconv)
% returns MI and distance dissimilarity matrices
% D = 1 - I(x;y) / H(x,y)

if ~islogical(deconv)
    deconv=logical(deconv);
end

p=sum(deconv)./size(deconv,1);
H=-(p.*log2(p)+(1-p).*log2(1-p));

p_joint=zeros(size(deconv,2),size(deconv,2),4);

p_joint(:,:,1)=double(deconv')*double(deconv);
p_joint(:,:,2)=double(~deconv')*double(deconv);
p_joint(:,:,3)=double(deconv')*double(~deconv);
p_joint(:,:,4)=double(~deconv')*double(~deconv);
p_joint=p_joint./sum(p_joint,3);

H_joint=-sum(p_joint.*log2(p_joint),3,'omitnan');
I=-H_joint+H+H';
I(isnan(I))=0;
D=1-I./H_joint;

I(diag(true(1,length(H))))=0;
D(diag(true(1,length(H))))=0;