function se=sem(x)
% standard error of the mean

se=std(x)./sqrt(size(x,1));