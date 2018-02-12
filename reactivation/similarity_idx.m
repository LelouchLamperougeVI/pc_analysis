function [M,mu,sd]=similarity_idx(M)

L=arrayfun(@(x) norm(M(x,:)),1:size(M,1));
L=L'*L;
M=M*M'./L;

iterations=100;
n=size(M,1)-1;
n=n*(n+1)/2;
coef=zeros(1,n*iterations);
for i=1:iterations
    shift=randi(size(M,1),1,size(M,2));
    temp=mat_circshift(M,shift);
    L=arrayfun(@(x) norm(temp(x,:)),1:size(temp,1));
    L=L'*L;
    temp=temp*temp'./L;
    temp=triu(temp,1);
    coef((i-1)*n+1:i*n)=temp(temp>0);
end

mu=mean(coef);
sd=std(coef);