function [Rsq,Rsq_adj,p]=lregress(x,y,n,h)
% Perform simple linear regression on variables x and y fitted to n-th
% order polynomial

if nargin < 3
    n=1;
end
if nargin < 4
    figure;
    h = gca;
end

p=polyfit(x,y,n);
yfit=polyval(p,x);
res=y-yfit;
SSres=sum(res.^2);
SStot=(length(y)-1)*var(y);

Rsq=1-SSres/SStot;

Rsq_adj=1 - SSres/SStot * (length(y)-1)/(length(y)-length(p));

plot(h,x,y,'k.');
hold on;
plot(h,x,yfit,'k');
axis square
title(h, ['R^{2}_{adj} = ' num2str(Rsq_adj)]);