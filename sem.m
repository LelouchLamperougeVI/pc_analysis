function se=sem(x)
% standard error of the mean

if size(x,1)==1
    x=x';
end

se=std(x)./sqrt(size(x,1));