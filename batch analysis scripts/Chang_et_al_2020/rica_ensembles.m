function rmd = rica_ensembles(X, clusts, varargin)
% Identify ensembles using reconstruction ICA
% Inputs:
%   X:
%   clusts:
% Name-Value pairs:
%   'norm':
%   'pca':
%   ''

ops = parse_ops(varargin);
tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro'); % measure explained variance

X(any(isnan(X), 2), :) = []; % get rid of movement epochs
if ops.norm
    X = (X - mean(X)) ./ std(X); % z-score X to null mean and unit variance
end

init = zeros(size(X, 2), length(clusts)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(clusts)
    init(clusts{ii}, ii) = 1;
end
init = init ./ vecnorm(init);

if ops.pca
    [coeff, score, latent] = pca(X);
    [~, l] = marchenkoPastur(X);
    sig_comps = latent > max(l);
    rmd.pca_ev = tev(X, coeff(:, sig_comps));
    X = score(:, sig_comps) * coeff(:, sig_comps)';
end

% md = rica(X, size(init, 2), 'initialtransformweights', init, 'contrastfcn', 'logcosh');
md = rica(X, sum(sig_comps));
W = md.TransformWeights;
W(:, sum(W) < 0) = -W(:, sum(W) < 0);
[~, idx] = max(W' * init);
W = W(:, idx);

rmd.sim = diag(W' * init);
rmd.ev = arrayfun(@(x) tev(X, W(:, x)), 1:size(W, 2));
rmd.tot_ev = tev(X, W);
rmd.z = X * W;
rmd.reconst = X * (W * W');
rmd.ops = ops;
rmd.W = W;
rmd.init = init;
rmd.X = X;

function ops = parse_ops(inputs)
ops.norm = true;
ops.pca = true;

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        case {'norm', 'normalize'}
            ops.norm = logical(inputs{count + 1});
        case {'pca', 'denoise'}
            ops.pca = logical(inputs{count + 1});
        otherwise
            error(['''' inputs{count} ''' is not a valid argument.']);
    end
    count = count + 2;
end