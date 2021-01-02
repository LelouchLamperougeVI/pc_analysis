function hiepi = hPICA(obj)
% Core hierarchical PCA/ICA (hPICA) method

tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro'); % measure explained variance

X = obj.twop.deconv;
if ~any(isnan(X), 'all')
    warning('Deconv does not contain NANs; possible omission of movement rejection.');
end

if isempty(obj.ensembles.clust)
    warning('No ensembles exist. Possible omission of ensemble detection.');
    hiepi = struct('pc', [], 'sim', [], 'ev', [], 'tot', [], 'W', [], 'init',...
            [], 'X', [], 'z', [], 'reconst', [], 'reactivations', [], 'psth', [],...
            'swr_react_strength', [], 'swr_mean_strength', [], 'swr_t', []);
    obj.hiepi = hiepi;
    return
end

% step 0. get rid of NANs and normalize to null mean, unitary variance
nan_idx = any(isnan(X), 2);
X(nan_idx, :) = []; % get rid of movement epochs
X = (X - mean(X)) ./ std(X); % z-score X to null mean and unit variance

% step 1. find initial estimate of basis vectors
init = zeros(size(X, 2), length(obj.ensembles.clust)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(obj.ensembles.clust)
    init(obj.ensembles.clust{ii}, ii) = 1;
end
init = init ./ vecnorm(init);

% step 2. estimate fraction of significant variance that need to be
% accounted for using Marchenko-Pastur theorem
[~, L] = eig(cov(X));
L = diag(L);
thres = march_pastur_lambda(X);
r_thres = size(X, 2) - sum(L(L < thres)); % the fraction of significant variance

% step 3. find eigenvectors of residual covariance matrix, get rid of the
% zero eigenvalue associated vectors and concatenate to basis estimate from
% hierarchical clustering to account for significant variance
resid = X - X * (init * init'); % residuals
[S, L] = eig(cov(resid));
[~, idx] = sort(diag(L));
S(:, idx(1:size(init, 2))) = [];
pc_comps = cat(2, init, S);
r = var(X * pc_comps);
[r((size(init, 2) + 1):end), idx] = sort(r((size(init, 2) + 1):end), 'descend');
pc_comps(:, (size(init, 2) + 1):end) = pc_comps(:, idx + size(init, 2));
r = cumsum(r);
W_ = pc_comps(:, r < r_thres);
if size(W_, 2) < size(init, 2)
    W_ = init;
end

X_ = X * W_;
md = rica(X_, size(W_, 2), 'InitialTransformWeights', eye(size(W_, 2)), 'ContrastFcn', 'logcosh');
W = md.TransformWeights' * W_';
W = W';
pc = W(:, 1:size(init, 2));

hiepi.pc = pc;
hiepi.sim = diag(pc' * init);
hiepi.ev = arrayfun(@(x) tev(X, W(:, x)), 1:size(W, 2));
hiepi.tot_ev = tev(X, W);
hiepi.W = W;
hiepi.init = init;

hiepi.X = nan(size(obj.twop.deconv));
hiepi.X(~nan_idx, :) = X;
hiepi.z = nan(size(obj.twop.deconv, 1), size(pc, 2));
hiepi.z(~nan_idx, :) = X * pc;
hiepi.reconst = nan(size(obj.twop.deconv));
hiepi.reconst(~nan_idx, :) = X * (W * W');

obj.hiepi = hiepi;
get_psth(obj);
detect_onsets(obj);

if ~isempty(obj.lfp.swr.swr)
    swr_response(obj);
    react_lfp(obj);
end


function l = march_pastur_lambda(X)
% Upper bound of eigenvalue defined by Marchenko-Pastur Law
dim = size(X);
l = (1 + sqrt(dim(2) / dim(1))) ^ 2;


function detect_onsets(obj)
% detect onset/offset of reactivations
z = movmedian(obj.hiepi.z, obj.twop.fs * obj.ops.hPICA.k);
obj.hiepi.md_filt_z = z;
thres = mean(z, 'omitnan') + obj.ops.hPICA.thres * std(z, 'omitnan');
ts = z > thres;

onset = z > (mean(z, 'omitnan') + obj.ops.hPICA.on_thres * obj.ops.hPICA.thres * std(z, 'omitnan'));

for e = 1:size(z, 2)
    heads = find(get_head(onset(:, e)));
    tails = get_head(onset(end:-1:1, e));
    tails = find(tails(end:-1:1));
    for ii = 1:length(heads)
        if ~any(ts(heads(ii):tails(ii), e))
            onset(heads(ii):tails(ii), e) = false;
        end
    end
    ts(:, e) = onset(:, e);
    
    gaps = obj.ops.hPICA.gap * obj.twop.fs;
    ts(:, e) = fill_gaps(ts(:, e), gaps);
    
    heads=find(get_head(ts(:, e)));
    tails=get_head(ts(end:-1:1, e));
    tails=find(tails(end:-1:1));

    if heads(1)==1
        heads(1)=[];
        tails(1)=[];
    end
    if tails(end)==length(ts)
        heads(end)=[];
        tails(end)=[];
    end

    obj.hiepi.reactivations{e}.peaks = zeros(length(heads),1);
    obj.hiepi.reactivations{e}.on = zeros(length(heads),1);
    obj.hiepi.reactivations{e}.dur = zeros(length(heads),1);
    for i = 1:length(heads)
        obj.hiepi.reactivations{e}.on(i) = obj.twop.ts(heads(i));
        [~, pks] = max(obj.hiepi.z(heads(i):tails(i), e));
        obj.hiepi.reactivations{e}.peaks(i) = obj.twop.ts(heads(i)+pks-1);
        obj.hiepi.reactivations{e}.dur(i) = obj.twop.ts(tails(i) + 1) - obj.twop.ts(heads(i));
    end
end


function get_psth(obj)
beh_dec = obj.analysis.original_deconv;
beh_dec = (beh_dec - min(beh_dec)) ./ range(beh_dec);
score = beh_dec * obj.hiepi.pc;
trial_bins = discretize(1:size(score,1), obj.analysis.behavior.trials_ts);
trial_bins(isnan(trial_bins)) = 0;
pos_bins = discretize(obj.analysis.behavior.unit_pos, length(obj.analysis.Pi));

psth = arrayfun(@(x) accumarray([pos_bins' trial_bins' + 1], score(:,x), [length(unique(pos_bins)) length(unique(trial_bins))], @mean), 1:size(score,2), 'uniformoutput', false);
psth = cellfun(@(x) x(:, 2:end), psth, 'uniformoutput', false);
obj.hiepi.psth = psth;


function swr_response(obj)
z = (obj.hiepi.z - mean(obj.hiepi.z, 'omitnan')) ./ std(obj.hiepi.z, 'omitnan');
one_sec = round(obj.twop.fs) * obj.ops.hPICA.swr_wdw;
t = knnsearch(obj.twop.ts, obj.lfp.swr.swr_on);
t = repmat(t, 1, one_sec * 2 + 1);
t = t + (-one_sec:one_sec);
t(any(t<1, 2), :) = [];
t(any(t>size(obj.twop.deconv, 1), 2), :) = [];

obj.hiepi.swr_react_strength = arrayfun(@(x) reshape(z(t, x), size(t)), 1:size(z, 2), 'UniformOutput', false);
obj.hiepi.swr_mean_strength = cellfun(@(x) mean(x, 'omitnan'), obj.hiepi.swr_react_strength, 'UniformOutput', false);
obj.hiepi.swr_mean_strength = cell2mat(obj.hiepi.swr_mean_strength')';

obj.hiepi.swr_t = linspace(-obj.ops.hPICA.swr_wdw, obj.ops.hPICA.swr_wdw, size(t, 2));


function react_lfp(obj)

for ii = 1:length(obj.hiepi.reactivations)
    trace = obj.lfp.swr.swr;
    p = envelope(trace, ceil(max(obj.lfp.fs ./ obj.lfp.ops.freqs('swr')) * 4), 'rms'); % half the nyquist frequency
    p = p .^ 2; % square for power
    p = (p - mean(p, 'omitnan')) ./ std(p, 'omitnan');
    one_sec = round(obj.lfp.fs) * obj.ops.hPICA.swr_wdw;
    t = knnsearch(obj.lfp.ts', obj.hiepi.reactivations{ii}.on);
    t = repmat(t, 1, one_sec * 2 + 1);
    t = t + (-one_sec:one_sec);
    t(any(t<1, 2), :) = [];
    t(any(t>length(trace), 2), :) = [];
    p = reshape(p(t), size(t));
    obj.hiepi.reactivations{ii}.swr_pw = p;
    
%     trace = obj.lfp.theta;
%     p = envelope(trace, ceil(max(obj.lfp.fs ./ obj.lfp.ops.freqs('theta')) * 4), 'rms'); % half the nyquist frequency
%     p = (p - mean(p, 'omitnan')) ./ std(p, 'omitnan');
%     one_sec = round(obj.lfp.fs) * obj.ops.hPICA.swr_wdw;
%     t = knnsearch(obj.lfp.ts', obj.hiepi.reactivations{ii}.on);
%     t = repmat(t, 1, one_sec * 2 + 1);
%     t = t + (-one_sec:one_sec);
%     t(any(t<1, 2), :) = [];
%     t(any(t>length(trace), 2), :) = [];
%     p = reshape(p(t), size(t));
%     obj.hiepi.reactivations{ii}.theta_pw = p;
%     
%     trace = obj.lfp.delta;
%     p = envelope(trace, ceil(max(obj.lfp.fs ./ obj.lfp.ops.freqs('delta')) * 4), 'rms'); % half the nyquist frequency
%     p = (p - mean(p, 'omitnan')) ./ std(p, 'omitnan');
%     one_sec = round(obj.lfp.fs) * obj.ops.hPICA.swr_wdw;
%     t = knnsearch(obj.lfp.ts', obj.hiepi.reactivations{ii}.on);
%     t = repmat(t, 1, one_sec * 2 + 1);
%     t = t + (-one_sec:one_sec);
%     t(any(t<1, 2), :) = [];
%     t(any(t>length(trace), 2), :) = [];
%     p = reshape(p(t), size(t));
%     obj.hiepi.reactivations{ii}.delta_pw = p;
%     
%     trace = obj.lfp.lgamma;
%     p = envelope(trace, ceil(max(obj.lfp.fs ./ obj.lfp.ops.freqs('lgamma')) * 4), 'rms'); % half the nyquist frequency
%     p = (p - mean(p, 'omitnan')) ./ std(p, 'omitnan');
%     one_sec = round(obj.lfp.fs) * obj.ops.hPICA.swr_wdw;
%     t = knnsearch(obj.lfp.ts', obj.hiepi.reactivations{ii}.on);
%     t = repmat(t, 1, one_sec * 2 + 1);
%     t = t + (-one_sec:one_sec);
%     t(any(t<1, 2), :) = [];
%     t(any(t>length(trace), 2), :) = [];
%     p = reshape(p(t), size(t));
%     obj.hiepi.reactivations{ii}.lgamma_pw = p;
%     
%     trace = obj.lfp.hgamma;
%     p = envelope(trace, ceil(max(obj.lfp.fs ./ obj.lfp.ops.freqs('hgamma')) * 4), 'rms'); % half the nyquist frequency
%     p = (p - mean(p, 'omitnan')) ./ std(p, 'omitnan');
%     one_sec = round(obj.lfp.fs) * obj.ops.hPICA.swr_wdw;
%     t = knnsearch(obj.lfp.ts', obj.hiepi.reactivations{ii}.on);
%     t = repmat(t, 1, one_sec * 2 + 1);
%     t = t + (-one_sec:one_sec);
%     t(any(t<1, 2), :) = [];
%     t(any(t>length(trace), 2), :) = [];
%     p = reshape(p(t), size(t));
%     obj.hiepi.reactivations{ii}.hgamma_pw = p;
end