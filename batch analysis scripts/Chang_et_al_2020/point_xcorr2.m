function [r, lags, p] = point_xcorr2(s1, s2, nans, ts, sig, wdw)

fs = 1 / median(diff(ts));
if nargin > 5
    wdw = wdw * fs;
else
    wdw = length(nans);
end

s1 = knnsearch(ts(:), s1(:));
s2 = knnsearch(ts(:), s2(:));

delay = s1 - s2';
edges = -(wdw - .5):(wdw - .5);
lags = -(wdw - 1):(wdw - 1);
r = histcounts(delay(:), edges);
r = fast_smooth(r, sig * fs, 2);

l1 = length(s1) / sum(~nans);
l2 = length(s2) / sum(~nans);

perms = 1e4;
stream = RandStream('dsfmt19937');
temp = rand(stream, [perms, sum(~nans), 2]) < permute([l1, l2], [1 3 2]);
null = false(perms, length(nans), 2);
null(:, ~nans, :) = temp;

r_null = zeros(perms, length(r));
for ii = 1:perms
    s1 = find(null(ii, :, 1));
    s2 = find(null(ii, :, 2));
    delay = s1(:) - s2(:)';
    
    r_null(ii, :) = histcounts(delay(:), edges);
end
r_null = fast_smooth(r_null, sig * fs, 2);

p = sum(r < r_null) ./ perms;