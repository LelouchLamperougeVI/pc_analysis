function hiepi = hPICA(obj)
% Core hierarchical PCA/ICA (hPICA) method

tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro'); % measure explained variance

X = obj.twop.deconv;
if ~any(isnan(X), 'all')
    warning('Deconv does not contain NANs; possible omission of movement rejection.');
end
if isempty(obj.ensembles.clust)
    warning('No ensembles exist. Possible omission of ensemble detection.');
    return
end

% step 0. get rid of NANs and normalize to null mean, unitary variance
X(any(isnan(X), 2), :) = []; % get rid of movement epochs
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

X_ = X * W_;
md = rica(X_, size(W_, 2), 'InitialTransformWeights', eye(size(W_, 2)), 'ContrastFcn', 'logcosh');
W = md.TransformWeights' * W_';
W = W';
pc = W(:, 1:size(init, 2));

hiepi.pc = pc;
hiepi.sim = diag(pc' * init);
hiepi.ev = arrayfun(@(x) tev(X, W(:, x)), 1:size(W, 2));
hiepi.tot_ev = tev(X, W);
hiepi.z = X * pc;
hiepi.reconst = X * (W * W');
hiepi.W = W;
hiepi.init = init;
hiepi.X = X;

obj.hiepi = hiepi;
get_psth(obj);


function l = march_pastur_lambda(X)
% Upper bound of eigenvalue defined by Marchenko-Pastur Law
dim = size(X);
l = (1 + sqrt(dim(2) / dim(1))) ^ 2;


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