function varargout=MPdist(M,N,lamb)
% Generate Marchenko Pastur distribution to threshold significance of
% eigenvalues for N variables and M observations
% Returns theoretical upper threshold and empirical probability for each eigenvalue in vector lamb

q=M/N;
if q<=1
    warning('Threshold prediction inaccurate due to q > 1 condition not satisfied');
end
a=(1+sqrt(1/q))^2+N^(-2/3);
b=(1-sqrt(1/q))^2;

varargout{1}=[b a];

if nargin<3 || nargout<2
    return;
end

pdf=@(q,a,b,lamb) q.*sqrt((a-lamb).*(lamb-b))./(2*pi.*lamb);
varargout{2}=pdf(q,a,b,lamb);