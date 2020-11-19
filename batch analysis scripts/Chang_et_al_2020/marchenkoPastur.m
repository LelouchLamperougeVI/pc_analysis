function [P, l] = marchenkoPastur(X, bins)
% Marchenko-Pastur distribution of random data

if nargin < 2
    bins = 100;
end

gamma = size(X,2) / size(X,1);
eig_range = [(1 + sqrt(gamma))^2 (1 - sqrt(gamma))^2];

l = linspace(eig_range(2), eig_range(1), bins);
P = sqrt( (eig_range(1) - l) .* (l - eig_range(2)) ) ./ (2 * pi * gamma .* l);