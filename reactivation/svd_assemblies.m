function assemblies=svd_assemblies(M,mu,sd,thres)

if nargin<4
    thres=2;
end
M=(M-mu)./sd>thres;
M=double(M);

[U,S,V]=svd(M);
s=diag(S);

sig=1;
while(s(sig)-s(sig+1)>2*(s(sig+1)-s(sig+2)))
    sig=sig+1;
end

assemblies=cell(1,sig);
for i=1:sig
    temp=zeros(size(S));
    temp(i)=s(i);
    temp=U*temp*V';
    assemblies{i}=logical(sum(temp));
end
