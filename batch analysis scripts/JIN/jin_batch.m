function res = jin_batch(root)

nperms = 1e4; % number of permutations for identifying speed cells
mvt_wdw = .5; % averaging window for movement epochs

sessions = dir(root);
idx = arrayfun(@(x) ~isempty(regexp(x.name, '\d', 'once')), sessions);
sessions = cat(1, sessions(idx).name);

if length(sessions) > 4
    error('Too many sessions, clean up the data.');
end
if ~all(ismember(sessions, '1234'))
    error('sessions must be either 1-through-3 or 1-through-4');
end

d_r1 = load(fullfile(root, '1', 'Plane1/deconv.mat'));
d_r1 = d_r1.deconv;
v_r1 = load(fullfile(root, '1', 'behavior.mat'));
v_r1 = v_r1.behavior.speed_raw;
if length(sessions) == 4
    d1 = load(fullfile(root, '2', 'Plane1/deconv.mat'));
    d2 = load(fullfile(root, '3', 'Plane1/deconv.mat'));
    s1 = load(fullfile(root, '2', 'behavior.mat'));
    s2 = load(fullfile(root, '3', 'behavior.mat'));
    deconv = cat(1, d1.deconv, d2.deconv);
    speed = cat(1, s1.behavior.speed_raw(:), s2.behavior.speed_raw(:));
    d_r2 = load(fullfile(root, '4', 'Plane1/deconv.mat'));
    d_r2 = d_r2.deconv;
    v_r2 = load(fullfile(root, '4', 'behavior.mat'));
    v_r2 = v_r2.behavior.speed_raw;
else
    deconv = load(fullfile(root, '2', 'Plane1/deconv.mat'));
    deconv = deconv.deconv;
    speed = load(fullfile(root, '2', 'behavior.mat'));
    speed = speed.behavior.speed_raw;
    d_r2 = load(fullfile(root, '3', 'Plane1/deconv.mat'));
    d_r2 = d_r2.deconv;
    v_r2 = load(fullfile(root, '3', 'behavior.mat'));
    v_r2 = v_r2.behavior.speed_raw;
end
tt = load(fullfile(root, '2', 'Plane1/timecourses.mat'));
fs = 1 / median(diff(tt.tcs.tt));

if mean(~sum(d_r1)) > .1
    warning('Rest1 data may be corrupt.')
end
if mean(~sum(deconv)) > .1
    warning('Run data may be corrupt.')
end
if mean(~sum(d_r2)) > .1
    warning('Rest2 data may be corrupt.')
end

d_r1 = ca_filt(d_r1);
d_r2 = ca_filt(d_r2);
deconv = ca_filt(deconv);

v_r1 = movsum(v_r1, round(fs * mvt_wdw)); % detect and reject moving epochs in rest
v_r2 = movsum(v_r2, round(fs * mvt_wdw));
d_r1(~~v_r1, :) = nan;
d_r2(~~v_r2, :) = nan;

d = deconv;
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

mua = cat(2, d_r1, deconv, d_r2) > 0;
mua = fast_smooth(mua, fs/2);
mua = reshape(mua, [size(deconv, 1), size(deconv, 2), 3]);
mua = cat(2, mean(mua(:, isspeed == -1, :), 2),...
    mean(mua(:, isspeed == 0, :), 2),...
    mean(mua(:, isspeed == 1, :), 2));
rho_mua = arrayfun(@(ii) corr(mua(:, :, ii)), 1:size(mua, 3), 'UniformOutput', false);
rho_mua = cat(3, rho_mua{:});

mu_fr = [mean(d_r1 > 0, 'omitmissing');...
    mean(deconv > 0, 'omitmissing');...
    mean(d_r2 > 0, 'omitmissing')]' .* fs;

res.isspeed = isspeed;
res.mu_fr = mu_fr;
res.rho_mua = rho_mua;