function md = sam_mape(x, n, trials, varargin)
% Maximum a posteriori estimation (aka Bayesian Decoding by
% neuroscientists) for Samsoon Inayat et al.
%
% Decoding error computed by Leave-one-out cross-validation (LOO-CV)
%
% Inputs:
%   x           distance/position vector
%   n           deconvolved firing rates matrix (r:time x c:neurons)
%   trials      trials vector of same length as x
% 'Name', 'Value' pairs:
%   'dt', 4             delta-t parameter for decoder (see Zhang et al. 1998)
%   'bins', 50          number of distance bins
%   'circ', false       compute circular error
%   'cv', true          leave-one-out cross-validation
%   'mle', false        use Maximum likelihood estimation instead (no prior)
%   'penalty', 1e-2     constant penalty for 0 firing rate in place of
%                       log(0) = -inf
%   'custom', vector    cross-validate on user specified indices
%
% Outputs:
%   md.decoded      decoded distance/position vector (cross-validated)
%   md.x            binned distance
%   md.err          mean absolute error as a function of position
%   md.oerr         overall mean error
%   md.mse          mean squared error (for debugging purposes)
%   md.ops          options used in getting this result
%
% Copyright (C) 2021 - HaoRan Chang
% SPDX-License-Identifier: GPL-2.0-or-later

ops = parse_inputs(varargin);

% remove NaNs
nan_idx = isnan(x);
x(nan_idx) = [];
trials(nan_idx) = [];
n(nan_idx, :) = [];
if ~isempty(ops.custom)
    ops.custom(nan_idx) = [];
end

% "straight up" the vectors
x = x(:);
trials = trials(:);

% there shall be no negative firing rate
n(n < 0) = 0;
n(isnan(n)) = 0;

% get rid of silent neurons
n(:, sum(n, 1) == 0) = [];

% bin distance
bin2cm = range(x) / ops.bins;
xbins = linspace(min(x), max(x), ops.bins + 1);
x = discretize(x, xbins);

decoded = zeros(size(x));
if ~isempty(ops.custom)
    decoded = decode(x(~ops.custom), n(~ops.custom, :), n(ops.custom, :), ops);
    x = x(ops.custom);
elseif ops.cv
    for k = unique(trials')
        loo = trials == k;
        decoded(loo) = decode(x(~loo), n(~loo, :), n(loo, :), ops);
    end
else
    decoded = decode(x, n, n, ops);
end

decoded = decoded(:);
if ops.circ
    err = min([mod(x - decoded, range(x)), mod(decoded - x, range(x))], [], 2);
else
    err = abs(x - decoded);
end
mse = mean((err .* bin2cm) .^ 2);
oerr = mean(err, 'omitnan');
err = accumarray(x, err, [ops.bins, 1], @mean);

cm = accumarray([x(:), decoded(:)], 1, [ops.bins, ops.bins]);
cm = cm ./ sum(cm, 2);

md.decoded = decoded .* bin2cm;
md.x = x .* bin2cm;
md.err = err .* bin2cm;
md.cm = cm;
md.oerr = oerr * bin2cm;
md.xbins = xbins;
md.mse = mse;
md.ops = ops;


function decoded = decode(x, train, test, ops)
% Gaussian smooth firing rate estimates
% Note: for 2p data, Gaussian smoothing greatly outperforms average
% windowing for whatever reason...

% make "stack"/lambda firing rates
Px = accumarray(x, 1, [ops.bins, 1]);

lambda = arrayfun(@(ii) accumarray(x, train(:, ii), [ops.bins, 1]), 1:size(train, 2), 'UniformOutput', false);
lambda = cell2mat(lambda);
lambda = lambda ./ Px;

Px = Px ./ length(x);

% summing window
kernel = ones(ops.dt, 1);
test = arrayfun(@(x) conv(test(:, x), kernel, 'same'), 1:size(test, 2), 'UniformOutput', false);
test = cell2mat(test);

% conduct decoding
n = permute(test, [3 1 2]); % dimensions: distance x time x neuron

lambda = fast_smooth(lambda, ops.sig);
lambda = lambda .* (1 - ops.penalty ./ range(lambda)) + ops.penalty;
lambda(lambda == 0) = nan;
lambda = permute(lambda, [1 3 2]);

ll = sum(n .* log(lambda) - ops.dt .* lambda, 3, 'omitnan');
if ~ops.mle
    ll = ll + log(Px);
end
[~, decoded] = max(ll, [], 1);


function ops = parse_inputs(inputs)
ops.bins = 50;
ops.dt = 4;
ops.sig = 5;
ops.circ = false;
ops.cv = true;
ops.mle = false;
ops.penalty = 1e-3;
ops.custom = [];

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        case {'bins', 'nbins', 'num_bins'}
            ops.bins = inputs{count+1};
        case {'dt', 'deltat', 'delta-t'}
            ops.dt = inputs{count+1};
        case {'sig', 'sd', 'smooth'}
            ops.sig = inputs{count+1};
        case {'circ', 'circular', 'mod', 'modulo'}
            ops.circ = logical(inputs{count+1});
        case {'cv', 'validate', 'cross-validate'}
            ops.cv = logical(inputs{count+1});
        case 'mle'
            ops.mle = logical(inputs{count+1});
        case {'penalty'}
            ops.penalty = inputs{count+1};
        case {'custom', 'user'}
            ops.custom = logical(inputs{count+1});
        otherwise
            error(['''' inputs{count} ''' is undefined.']);
    end
    count = count + 2;
end


function smoothed=fast_smooth(data,sig,dim)
% Faster 1d gaussian kernel smoothing
% Accounts for edge underestimation and NaNs
% Usage: 
%   data: matrix of size n x m where rows contain observations and columns
%      contain variables (can also be vector)
%   sig: kernel standard dev
%   dim: dimension of observations (default 1)
%
%   smoothed: smoothed data
%
% by HaoRan Chang

if sig==0
    smoothed=data;
    return
end

if nargin==2
    dim=1;
end

switch dim
    case 1
    case 2
        data=data';
    otherwise
        error('WTF dude?')
end

dataSize=size(data,1);
kernelSize=ceil(10*sig);
kernelSize=kernelSize-~mod(kernelSize,2); %avoid even kernel sizes
alpha=(kernelSize-1)/sig/2;
kernel=gausswin(kernelSize,alpha);
kernel=kernel./sum(kernel);

taper=zeros(kernelSize,1);
isnan_idx = isnan(data); % now generalized
nan_idx = [repmat(taper,1,size(data,2)); ~isnan_idx; repmat(taper,1,size(data,2))];
nan_idx = nan_idx(:);
data=[repmat(taper,1,size(data,2));data;repmat(taper,1,size(data,2))];
data=data(:);
data(isnan(data))=0;

smoothed=conv(data,kernel,'same');
smoothed=reshape(smoothed,dataSize+kernelSize*2,[]);
smoothed = smoothed./ reshape( conv(nan_idx - 1, kernel,'same') + sum(kernel), dataSize+kernelSize*2, [] ); % new performance oriented (^old) with generalization
smoothed([1:kernelSize end-kernelSize+1:end],:)=[];
smoothed(isnan_idx) = nan;

switch dim
    case 1
    case 2
        smoothed=smoothed';
end