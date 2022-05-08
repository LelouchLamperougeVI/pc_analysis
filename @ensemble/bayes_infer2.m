function [decoded, P, pos, err] = bayes_infer2(obj,varargin)
% New Bayesian decoding routine
% Wrapper for sam_mape
%
% Usage
% 'name', value pairs with default option
%     'dt', 1             Decoding time window in seconds (see Zhang et al. 1998)
%     'smooth', 0.2       Gaussian smoothing factor in seconds
%     'bins', 50          Number of position bins
%     'circ', true        Compute circular distance for error calculation
%     'cv', 'loo'         Cross-validation method. By default, use leave-one-out cross-validation ('loo'). Alternatively, can use even-odd ('eo')
%     'mle', false        Run maximum likelihood estimation instead of maximum a posteriori
%     'penalty', eps      Constant penalty factor in case of log(0) = -inf (i.e. 0 firing rates)
%     
% outputs
%     decoded     Decoded position vector
%     P           Log-likelihood matrix of size bins X positions
%     pos         Actual animal position vector
%     err         Decoding error as a function of position with (2 columns matrix with [error, standard error of the mean])

if isempty(obj.analysis)
    error('run analysis first you dumb dumb');
end

ops = parse_inputs(varargin);
analysis = obj.analysis;

x = analysis.behavior.unit_pos(:);
n = analysis.original_deconv;
trials = analysis.behavior.trials_ts(:);
trials = discretize(1:length(x), trials);
trials = trials(:);

if ops.rm_mvt
    thres = noRun(analysis.behavior.unit_vel);
    thres=( abs(analysis.behavior.unit_vel) > thres ) &...
        (analysis.behavior.trials(1) < analysis.behavior.frame_ts &...
        analysis.behavior.trials(end) > analysis.behavior.frame_ts);
    
    x = x(thres);
    n = n(thres, :);
    trials = trials(thres);
end

if strcmpi(ops.cv, 'loo')
    ops.cv = [];
elseif strcmpi(ops.cv, 'eo')
    ops.cv = ~~mod(trials, 2);
else
    error('invalid cross-validation method');
end

ops.dt = round(analysis.fs * ops.dt);
ops.smooth = round(analysis.fs * ops.smooth);

md = sam_mape(x, n, trials, 'dt', ops.dt, 'smooth', ops.smooth, 'bins', ops.bins, ...
    'circ', ops.circ, 'cv', true, 'mle', ops.mle, 'penalty', ops.penalty, 'custom', ops.cv);

decoded = md.decoded;
P = md.ll;
pos = md.x;
err = [md.err(:), md.sem(:)];

obj.bayes = md;


function ops = parse_inputs(inputs)
ops.dt = 1; % all in seconds
ops.smooth = .2;
ops.bins = 50;
ops.circ = true;
ops.cv = 'loo';
ops.mle = false;
ops.penalty = eps;

ops.rm_mvt = true;

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        case {'dt'}
            ops.dt = inputs{count + 1};
        case {'smooth', 'sd'}
            ops.smooth = inputs{count + 1};
        case {'bins', 'bin'}
            ops.bins = inputs{count + 1};
        case {'circ', 'circular'}
            ops.circ = inputs{count + 1};
        case {'cv', 'cross validate', 'cross-validate'}
            ops.cv= inputs{count + 1};
        case {'mle', 'maximum likelihood', 'maximum likelihood estimation'}
            ops.mle = inputs{count + 1};
        case {'penalty'}
            ops.penalty = inputs{count + 1};
        case {'mvt', 'rm_mvt'}
            ops.rm_mvt = inputs{count + 1};
        otherwise
            error('invalid input');
    end
    
    count = count + 2;
end