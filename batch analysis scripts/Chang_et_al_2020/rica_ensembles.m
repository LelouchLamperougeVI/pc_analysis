function rmd = rica_ensembles(X, clusts)
% Identify ensembles using reconstruction ICA
% Inputs:
%   X:          t x n data matrix
%   clusts:     cell vector of ensemble clusters

tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro'); % measure explained variance

% step 0. get rid of NANs and normalize to null mean, unitary variance
X(any(isnan(X), 2), :) = []; % get rid of movement epochs
X = (X - mean(X)) ./ std(X); % z-score X to null mean and unit variance

% step 1. find initial estimate of basis vectors
init = zeros(size(X, 2), length(clusts)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(clusts)
    init(clusts{ii}, ii) = 1;
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

rmd.pc = pc;
rmd.sim = diag(pc' * init);
rmd.ev = arrayfun(@(x) tev(X, W(:, x)), 1:size(W, 2));
rmd.tot_ev = tev(X, W);
rmd.z = X * pc;
rmd.reconst = X * (W * W');
rmd.W = W;
rmd.init = init;
rmd.X = X;


function l = march_pastur_lambda(X)
% Upper bound of eigenvalue defined by Marchenko-Pastur Law
dim = size(X);
l = (1 + sqrt(dim(2) / dim(1))) ^ 2;