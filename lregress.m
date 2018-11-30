function [Rsq,Rsq_adj,p]=lregress(x,y,n)
% Perform simple linear regression on variables x and y fitted to n-th
% order polynomial

if nargin < 3
    n=1;
end

p=polyfit(x,y,n);
yfit=polyval(p,x);
res=y-yfit;
SSres=sum(res.^2);
SStot=(length(y)-1)*var(y);

Rsq=1-SSres/SStot;

Rsq_adj=1 - SSres/SStot * (length(y)-1)/(length(y)-length(p));

figure;hold on;
plot(x,y,'.');
plot(x,yfit);
title(['R^{2}_{adj} = ' num2str(Rsq_adj)]);