function [p, chi] = chisq2(A, B)
% Conduct a two-sample chi squared test
% A and B have a finite set of unique identifiers
% By default, chi squared test is two-tailed

classes = unique(cat(1, A(:), B(:)));
A = categorical(A(:), classes);
B = categorical(B(:), classes);
classes = categorical(classes);

o = cat(1, histcounts(A, classes), histcounts(B, classes));
e = sum(o, 1) ./ sum(o, 'all') .* sum(o, 2);

chi = sum((o - e).^2 ./ e, 'all');
df = prod(size(o) - 1);

p = 1 - chi2cdf(chi, df);