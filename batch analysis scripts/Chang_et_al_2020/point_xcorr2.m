function [r, lags, p] = point_xcorr2(s1, s2, nans, ts, sig, wdw)

fs = 1 / median(diff(ts));
if nargin > 5
    wdw = wdw * fs;
else
    wdw = length(nans);
end

s1 = knnsearch(ts(:), s1(:));
s2 = knnsearch(ts(:), s2(:));

N = length(ts);
delay = s1 - s2';
edges = -(N - .5):(N - .5);
lags = -(N - 1):(N - 1);
% r = histcounts(delay(:), edges) ./ numel(delay);
r = histcounts(delay(:), edges);
r = fast_smooth(r, sig * fs, 2);

l1 = length(s1) / sum(~nans);
l2 = length(s2) / sum(~nans);

perms = 1e3;
stream = RandStream('dsfmt19937');
temp = rand(stream, [perms, sum(~nans), 2]) < permute([l1, l2], [1 3 2]);
null = false(perms, N, 2);
null(:, ~nans, :) = temp;

r_null = zeros(perms, length(r));
for ii = 1:perms
    s1 = find(null(ii, :, 1));
    s2 = find(null(ii, :, 2));
    delay = s1(:) - s2(:)';
    
%     r_null(ii, :) = histcounts(delay(:), edges) ./ numel(delay);
    r_null(ii, :) = histcounts(delay(:), edges);
end
r_null = fast_smooth(r_null, sig * fs, 2);

p = sum(r < r_null) ./ perms;

r(abs(lags) > wdw) = [];
p(abs(lags) > wdw) = [];
lags(abs(lags) > wdw) = [];
lags = lags ./ fs;