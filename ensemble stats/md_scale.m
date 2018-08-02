function [X,V,e,strain]=md_scale(D,plotFlag)
% Given similarity/dissimilarity matrix D, perform classical 
% multidimensional scaling and return eigen vectors and values V and e
if nargin <2
    plotFlag=false;
end

centering_matrix=@(n) eye(n)-(1/n).*ones(n);

C=centering_matrix(length(D));

D=D.^2;
B=C*(-.5.*D)*C;

[V,e]=eig((B+B')./2); %guard against spurious eigenvalues? just don't question it...
[~,idx]=sort(diag(e),'descend');
V=V(:,idx);
e=e(idx,idx);

X=V*sqrt(e);
[~,idx]=max(abs(X));
X=X.*sign(X(idx+size(X,1).*(0:size(X,2)-1)));

strain=get_strain(B,V,e);
if plotFlag
    plot(strain);
    xlabel('number eigen components');
    ylabel('loss');
end

function strain=get_strain(B,V,e)

b_sqrt=B.^2;
b_sqrt=sum(b_sqrt(:));

strain=zeros(length(e),1);
for i=1:length(e)
    X=V(:,1:i)*sqrt(e(1:i,1:i));
    X=X*X';
    X=(B-X).^2;
    X=sum(X(:));
    strain(i)=sqrt(X/b_sqrt);
end