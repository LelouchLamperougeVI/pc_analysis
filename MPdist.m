function e=MPdist(M,N,lamb)
% Generate Marchenko Pastur distribution to threshold significance of
% eigenvalues for N variables and M observations
% Returns empirical probability for each eigenvalue in vector lamb

q=M/N;
a=(1+sqrt(1/q))^2;
b=(1-sqrt(1/q))^2;

pdf=@(q,a,b,lamb) q.*sqrt((a-lamb).*(lamb-b))./(2*pi.*lamb);
e=pdf(q,a,b,lamb);