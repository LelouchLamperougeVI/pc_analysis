function n=fd_bins(data)
% Apply the Freedman-Diaconis rule to determine the "optimal" number of
% bins

data(isinf(data))=nan;
n=sum(~isnan(data));
iq_range=iqr(data);
n=2.*iq_range.*(n.^(-1/3));

n=floor(range(data)./n);