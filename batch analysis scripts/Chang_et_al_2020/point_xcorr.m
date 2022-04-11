function [r, lags, p] = point_xcorr(s1, s2, N, fs)
% Model the crosscorrelogram between two point processes

s1 = s1(:); s2 = s2(:);
% estimation of true sampling rate given the approx. fs
% this stupid step is necessary because I am a fucking dumbass
isint = @(x) all(rem(x(:), 1) == 0);
if ~isint(s1) || ~isint(s2)
    if nargin < 4
        error('s1, s2 are timestamps, but no sampling rate was given');
    end
    
    temp1 = s1 - s1';
    temp1 = abs(triu(temp1, 1));
    temp1 = temp1(~~temp1);
    temp2 = s2 - s2';
    temp2 = abs(triu(temp2, 1));
    temp2 = temp2(~~temp2);
    temp = cat(1, temp1(:), temp2(:));
    
    temp = sort(temp);
    temp = temp(1:min([50, length(temp)]));
    fs = round(temp .* fs);
    fs = median(fs ./ temp);
    
    disp(['Sampling rate estimate: ' num2str(fs)]);
else
    fs = 1;
end

delay = s1 - s2';
lags = -(N - 1):(N - 1);
lags = lags ./ fs;
edges = -(N - .5):(N - .5);
edges = edges ./ fs;
r = histcounts(delay(:), edges);

l1 = length(s1) / N;
l2 = length(s2) / N;
pp = @(lamb1, lamb2, k, N) (lamb1 .* lamb2 .* N) .^ k .* exp(-lamb1 .* lamb2 .* N) ./ factorial(k);
p = pp(l1, l2, r, N - lags);