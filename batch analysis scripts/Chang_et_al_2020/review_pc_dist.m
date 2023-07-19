clear all
ee = load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'pc_list', 'whole_stack', 'pc_width');
rsc = load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'pc_list', 'whole_stack', 'pc_width');

whole_stack = cat(2, rsc.whole_stack, ee.whole_stack);
pc_list = cat(2, rsc.pc_list, ee.pc_list);
pc_width = cat(2, rsc.pc_width, ee.pc_width);
pc_width = cat(2, pc_width{:})';
pc_width(cellfun(@isempty, pc_width)) = [];

stack = arrayfun(@(ii) whole_stack{ii}(:, pc_list{ii}), 1:length(pc_list), 'UniformOutput', false);
stack = cat(2, stack{:});

[~, order] = max(stack);
[~, order] = sort(order);

figure
subplot(2, 2, [1 3]);
imagesc(stack(:, order)')
colormap viridis

loc = cat(1, pc_width{:});
loc = loc(:, 2);

[~, loc] = max(stack);

belt;
belt_idx = ~~conv(belt_idx, [1 1 1 1 1], 'same');

O = histcounts(loc);
x0 = O(belt_idx);
O = O(~belt_idx);

O(1) = []; O(end) = [];

n = range(loc) + 1;
% chi2 = sum((O - E).^2 ./ E);

chi = @(x) sum(([O(:); x(:)] - (sum(O) + sum(x)) / n).^2 ./ ((sum(O) + sum(x)) / n));
est = fmincon(chi, x0(:), eye(length(x0)), x0);

real = histcounts(loc);
estimated = real;
estimated(belt_idx) = est;
real = real ./ sum(real);
estimated = estimated ./ sum(estimated);

subplot(2, 2, 2)
plot(real)
hold on
plot(estimated)
ylim([0 .15])

frac = sum(x0(:) - est) / length(loc);
subplot(2, 2, 4);
pie([1 - frac, frac]);