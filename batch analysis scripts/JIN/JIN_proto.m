% load('/home/loulou/Documents/Data/uleth/rsc036/2017_08_03/2/Plane1/deconv.mat')
% load('/home/loulou/Documents/Data/uleth/rsc036/2017_08_03/2/behavior.mat')
% load('/home/loulou/Documents/Data/uleth/rsc036/2017_08_03/2/Plane1/timecourses.mat')
% 
% f = logspace(-2, 1, 2^8);
% fs = 1 / mean(diff(tcs.tt));
% d = ca_filt(deconv);
% d = mean(d > 0, 2);
% d = d .* fs; % normalise to mean #spikes per second
% d = fast_smooth(d, fs/2);
% 
% wdw = round(fs*200);
% nover = round(wdw * .75);
% p = pwelch(d, hann(wdw), nover, f, fs, 'psd');
% 
% figure
% semilogx(f, p)
% 
% edges = logspace(-3, 0, 51);
% v = abs(behavior.speed_raw);
% v(v < 1e-3) = 1e-3;
% p = accumarray([discretize(d, 200), discretize(v, edges)], 1);
% % p = accumarray([discretize(d, 100), discretize(v, 100)], 1);
% p = p ./ sum(p);
% p(isnan(p)) = 0;
% figure
% imagesc(imgaussfilt(p, 5))

clear all

root = '/home/loulou/Documents/Data/uleth/rsc036';

regions = dir(root);
regions = {regions(3:end).name};

results = {};
for r = regions
    sessions = dir(fullfile(root, r{1}));
    sessions = {sessions(3:end).name};

    for s = sessions
        try
            res = jin_batch(fullfile(root, r{1}, s{1}));
        catch
            continue
        end
        isspeed = repmat(res.isspeed(:) + 2, [3, 1]);
        block = repelem(1:3, size(res.mu_fr, 1));
        mu_fr = accumarray([isspeed(:), block(:)], res.mu_fr(:), [3, 3], @mean);
        idx = triu(ones(length(res.rho_mua)), 1);
        rho = [];
        for ii = 1:3
            rho_mua = res.rho_mua(:, :, ii);
            rho = [rho; rho_mua(~~idx)];
        end
        isspeed = histcounts(res.isspeed, -1:2) ./ length(res.isspeed);
        results = [results; {r{1}, mu_fr(:)', rho(:)', isspeed(:)'}];

        disp(fullfile(root, r{1}, s{1}))
    end
end

%%
isspeed = cat(1, results{:, 4});
x = results(:, 1);
x = repmat(x, [1, 3]);
y = repmat({'speed-off', 'none', 'speed-on'}, [size(x, 1), 1]);

g = gramm('x', x(:), 'y', isspeed(:), 'color', y(:));
g.stat_boxplot()
figure
g.draw()

mu_fr = cat(1, results{:, 2});
mu_fr = mu_fr(:, 4:6) ./ mu_fr(:, 1:3);
x = results(:, 1);
x = repmat(x, [1, 3]);
y = repmat({'speed-off', 'none', 'speed-on'}, [size(x, 1), 1]);

g = gramm('x', x(:), 'y', mu_fr(:), 'color', y(:));
g.stat_boxplot()
figure
g.draw()

mu_fr = cat(1, results{:, 2});
mu_fr = mu_fr(:, 4:6) ./ mu_fr(:, 7:9);
x = results(:, 1);
x = repmat(x, [1, 3]);
y = repmat({'speed-off', 'none', 'speed-on'}, [size(x, 1), 1]);

g = gramm('x', x(:), 'y', mu_fr(:), 'color', y(:));
g.stat_boxplot()
figure
g.draw()

rho_mua = cat(1, results{:, 3});
rho_mua = rho_mua(:, 1:3);
x = results(:, 1);
x = repmat(x, [1, 3]);
y = repmat({'-1 0', '-1 1', '0 1'}, [size(x, 1), 1]);

g = gramm('x', x(:), 'y', rho_mua(:), 'color', y(:));
g.stat_boxplot()
figure
g.draw()

rho_mua = cat(1, results{:, 3});
rho_mua = rho_mua(:, 7:9);
x = results(:, 1);
x = repmat(x, [1, 3]);
y = repmat({'-1 0', '-1 1', '0 1'}, [size(x, 1), 1]);

g = gramm('x', x(:), 'y', rho_mua(:), 'color', y(:));
g.stat_boxplot()
figure
g.draw()

%% identify speed modulated cells
nperms = 1e4;

d = ca_filt(deconv);
v = abs(speed);
d = single(zscore(d));
v = single(zscore(v));

rho = @(x, y) x' * y ./ length(x);
rsq = rho(v, d)'.^ 2;
rsq_n = zeros(size(deconv, 2), nperms);
for ii = 1:nperms
    rsq_n(:, ii) = rho(circshift(v, randi(size(deconv, 1))), d) .^ 2;
end
p = sum(rsq <= rsq_n, 2) ./ nperms;
p(p == 0) = 1 / (nperms + 1);

rsq = rho(v, d)';
isspeed = (p < .01) .* sign(rsq); % speed-on, nothing, speed-off

figure
boxplot(mean(deconv > 0, 'omitmissing') ./ mean(d_r1 > 0, 'omitmissing'), isspeed)

%%
[~, order] = sort(rsq);

figure
ax(1) = subplot(4, 1, 1:3);
imagesc(zscore(fast_smooth(deconv(:, order), 10))')
ax(2) = subplot(4, 1, 4);
plot(speed)
linkaxes(ax, 'x')

%%
x = repmat(discretize(abs(speed), 100), [1, size(deconv, 2)]);
y = repmat(1:size(deconv, 2), [size(deconv, 1), 1]);
stack = accumarray([x(:), y(:) ], deconv(:), [], @mean);
stack = (stack - min(stack)) ./ range(stack);
stack = zscore(stack);
stack = fast_smooth(stack, 5);

figure
[~, order] = sort(rsq);
imagesc(stack(:, order)')

%%
f = logspace(-2, 0, 2^8);
fs = 1 / mean(diff(tcs.tt));

mua = ca_filt(deconv) > 0;
mua = fast_smooth(mua, fs/2);
mua = [mean(mua(:, isspeed == 1), 2),...
    mean(mua(:, isspeed == -1), 2),...
    mean(mua(:, isspeed == 0), 2)];
mua = mua .* fs;

wdw = round(fs*200);
nover = round(wdw * .75);
psd = pwelch(mua(3e3:7e3,:), hann(wdw), nover, f, fs, 'psd');

figure
semilogx(f, psd)

%%
% 
% [~, order] = sort(corr(v, fast_smooth(deconv, 10)) .^ 2, 'descend');
% rsq = zeros(length(order), 1);
% for ii = 1:length(rsq)
%     md = fitlm(fast_smooth(deconv(:, order(1:ii)), 10), v);
%     rsq(ii) = md.Rsquared.Ordinary;
% end