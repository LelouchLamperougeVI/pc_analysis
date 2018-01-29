function C=simple_conv(A,k)
% 1D convolution that also accepts matrices
% Signal and kernel over columns
if size(A,1)<size(k,1)
    error('kernel longer than signal');
end

C=zeros(size(A,1)+size(k,1)-1,size(A,2));

A=[zeros(size(k,1)-1,size(k,2));A;zeros(size(k,1)-1,size(k,2))];
k=flipud(k);
n=size(k,1);
for i=1:size(C,1)
    C(i,:)=sum(A(i:n+i-1,:).*k);
end