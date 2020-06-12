function se=sem(x,dim)
% standard error of the mean

if nargin<2
    dim = 1;
end

se=std(x,0,dim,'omitnan')./sqrt(sum(~isnan(x),dim));