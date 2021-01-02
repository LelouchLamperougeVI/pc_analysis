function [r, lags, p] = bxcorr(X, Y, t)
% MATLAB built-in xcorr function kinda sucks
% Here, we compute the real CCF
if ~isvector(X) || ~isvector(Y)
    error('Inputs must be vectors');
end
if length(X) ~= length(Y)
    error('Inputs vectors must be of equal lengths');
end
X = X(:);
Y = Y(:);

N = length(X);
X_ = X - mean(X, 'omitnan');
Y_ = Y - mean(Y, 'omitnan');
lags = (-(N - 1) : (N - 1))';
k = abs(lags);
N_ = N - k;
Sx = sqrt(conv(X_.^2, ones(N, 1)) ./ N_);
Sy = sqrt(conv(Y_.^2, ones(N, 1)) ./ N_);
r = conv(X_, Y_(end:-1:1)) ./ N_ ./ Sx ./ Sy;

lags = lags .* median(diff(t));
p = 2 ./ sqrt(N);