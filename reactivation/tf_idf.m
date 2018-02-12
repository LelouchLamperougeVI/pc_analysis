function M=tf_idf(X)
% calcium events normalization using search algorithm described in Reid et
% al. 2016

tf=1./sum(X,2);
idf=log(size(X,1)./sum(X));
M=tf*idf.*X;

M(isnan(M))=0;