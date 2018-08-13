function [I,D]=get_mi_knn(deconv)

I = zeros(size(deconv,2));
ds=[1;1];
co=IKCCA_initialization(1);
% co=IShannon_HShannon_initialization(1);

for i=1:size(deconv,2)
    parfor j=i+1:size(deconv,2)
        I(i,j)=IKCCA_estimation(deconv(:,[i j])',ds,co);
%         I(i,j)=IShannon_HShannon_estimation(deconv(:,[i j])',ds,co);
    end
    I(:,i)=I(i,:);
    i
end

D=max(I(:))-I;

D=-D.*(1-diag(ones(1,size(deconv,2))));