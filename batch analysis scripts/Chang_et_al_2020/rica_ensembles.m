function [W, z, reconstructed, init] = rica_ensembles(X, clusts, varargin)
% Identify ensembles using reconstruction ICA
% Inputs:
%   X:
%   clusts:
% Name-Value pairs:
%   'norm':
%   'pca':
%   ''

ops = parse_ops(varargin);

tev = @(X, W) norm(X * (W * W'), 'fro') ./ norm(X, 'fro');

X(any(isnan(X), 2), :) = [];
normX = (X - mean(X)) ./ std(X); % z-score X to null mean and unit variance

init = zeros(size(X, 2), length(clusts)); % construct initial estimate of mixing matrix W tilde
for ii = 1:length(clusts)
    init(clusts{ii}, ii) = 1;
end
init = init ./ vecnorm(init);

X_hat = normX - normX * (init * init');
[coeff, ~, latent] = pca(X_hat);
[~, l] = marchenkoPastur(X);
coeff = coeff(:, latent <= max(l));
init = cat(2, init, coeff);
init(:, size(init, 1) + 1:end) = [];
while tev(normX, init(:, 1:end-1)) > ev_thres && size(init, 2) > length(clusts)
    init(:,end) = [];
end

md = rica(normX, size(init, 2), 'initialtransformweights', init, 'contrastfcn', 'logcosh');
W = md.TransformWeights(:, 1:length(clusts));

if nargout > 1
    z = normX * W;
    reconstructed = normX * (W * W');
end

function ops = parse_ops(inputs)
ops.norm = true;
ops.pca = true;