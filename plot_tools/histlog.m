function h=histlog(X,n)
% logscale histogram

if nargin<2
    n=50; %default
end

x=logspace(floor(log10(min(X))),ceil(log10(max(X))),n);
h=histogram(X,x);
set(gca,'xscale','log');