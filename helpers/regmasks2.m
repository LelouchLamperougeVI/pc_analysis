function [I, J, reg] = regmasks2(m1, m2, neighbourhood)
% A simplified version of regmasks omitting rotations (only translations).
% Does the following:
%   1. Register m2 (moving) against m1 (fixed) using the convolution
%      theorem trick.
%   2. Calculate the Jaccard index for each pair of ROIs between m1 and m2
%   3. Use k-means to separate the Jaccard indices into two clusters of
%      true and false positives
%
% Usage:
%   [I, J, reg] = regmasks2(m1, m2, neighbourhood)
%
% Inputs:
%   m1:             the reference mask
%   m2:             the mask to be registered (moving)
%   neighbourhood:  search radius (in fractions) - default = .05
%
% Outputs:
%   I:      Vector of same length as number of ROIs in m1, containing
%           the IDs of matching ROIs in m2. Neurons in m1 with no matching
%           ROI in m2 have values NAN.
%   J:      Vector of Jaccard distances for each ROI. Same length as I.
%   reg:    Registered m2 mask.

if nargin < 3
    neighbourhood = .05;
end

neighbourhood = min(round(size(m1) .* neighbourhood));

b1 = double(~~m1);
b2 = double(~~m2);
if numel(m1) ~= numel(m2)
    error('mismatching dimensions between masks');
end

grid = 2 .* size(b1) - 1;
x = (1:grid(2)) - size(b1, 2);
y = (1:grid(1)) - size(b1, 1);

a = padarray(b1, size(b2) - 1, 'post');
b = padarray(b2(end:-1:1, end:-1:1), size(b1) - 1, 'post');
c = ifft2( fft2(a) .* fft2(b) ); % ye olde convolution theorem

c = c(abs(y) <= neighbourhood, abs(x) <= neighbourhood);
x = x(abs(x) <= neighbourhood);
y = y(abs(y) <= neighbourhood);

[~, ind] = max(c(:));
[py, px] = ind2sub(size(c), ind);

reg = imtranslate(m2, [x(px), y(py)]);

target = sparse(m1(:)) == (1:range(m1(:))); % I'm lazy so sparse arrays
match = sparse(reg(:)) == (1:range(reg(:)));

J = target' * match;
J = J ./ ( sum(target)' + sum(match) - J );
[J, I] = max(J, [], 2);
J = full(J);

k = kmeans(J, 2);
if mean(J(k == 1)) > mean(J(k == 2))
    k = ~(k - 1) + 1;
end

J(k == 1) = nan;
I(k == 1) = nan;
