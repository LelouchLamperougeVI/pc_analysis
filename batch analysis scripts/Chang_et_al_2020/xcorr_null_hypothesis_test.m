% Given two i.i.d. random point processes and a sufficiently small delta-t,
% the probability of observing co-occurence of the two events is given by
% a Poisson distribution at a rate equivalent to the product of the rate of
% the individual point processes.
% We test the validity of this hypothesis here.
% Verdict: yup, that's indeed the case

l1 = .6;
l2 = .05;
N = .5e3;
boot = 1e4;

spk = rand(boot, N, 2) < permute([l1, l2], [1 3 2]);
co = spk(:, :, 1) & spk(:, :, 2);
co = sum(co, 2);

pp = @(lamb1, lamb2, k, N) (lamb1 * lamb2 * N) .^ k .* exp(-lamb1 * lamb2 * N) ./ factorial(k);

figure
histogram(co, 'Normalization', 'probability');
hold on
k = 0:100;
plot(k, pp(l1, l2, k, N));


%%
dt = .001;
l1 = 100; %1.2;
l2 = 140; %5.6;
T = 1;
wdw = .002;
perms = 1e4;

N = T / dt;

wdw = wdw / dt;
edges = -(N - wdw/2):wdw:(N - wdw/2);
lags = -(N - wdw):wdw:(N - wdw);

r = zeros(perms, length(lags));
for ii = 1:perms
    spk = rand(N, 2) < ([l1, l2] .* dt);
    s1 = find(spk(:, 1));
    s2 = find(spk(:, 2));

    delay = s1(:) - s2(:)';

    r(ii, :) = histcounts(delay(:), edges);
end

p_boot = arrayfun(@(x) histcounts(r(:, x), -0.5:400.5), 1:size(r, 2), 'UniformOutput', false);
p_boot = cell2mat(p_boot') ./ perms;

lambda = l1 * l2 * dt^2;
n = N - abs(lags(:));
n = (2 .* n - wdw) .* wdw ./ 2;
k = 0:400;
% p = k .* log(lambda .* n .* wdw) - lambda .* n .* wdw - gammaln(k + 1);
p = k .* log(lambda .* n) - lambda .* n - gammaln(k + 1);
p = exp(p);

figure
plot(p(lags == 0, :))
hold on
plot(p_boot(lags == 0, :))
% p = cumsum(p, 2);

% figure
% % imagesc('xdata', lags, 'ydata', k, 'cdata', p');
% imagesc('xdata', lags, 'ydata', k, 'cdata', p');
% hold on
% plot(lags, r, 'r');



%%
flat = rand(perms, N) < .01;
flat = sum(flat, 2);
binned = rand(perms, N/100, 100) < .01;
binned = sum(flat, [2, 3]);

