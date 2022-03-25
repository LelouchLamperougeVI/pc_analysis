function [decoded, P, pos, err] = bayes_paired(a1, a2, ROIs, varargin)
% sam_mape wrapper for cross days stack
% ROI index maps a1 onto a2
% i.e., a1.deconv(ROIs) -> a2.deconv
% Train with a1, test on a2
% Rest of it is the same as bayes_infer2...

if abs(a1.fs - a2.fs) > 1e-2
    error('The two sessions were acquired at different sampling rates...');
end

ops = parse_inputs(varargin);

cv_idx = repelem([false; true], [length(a1.behavior.unit_pos); length(a2.behavior.unit_pos)]);
x = cat(1, a1.behavior.unit_pos(:), a2.behavior.unit_pos(:));
n = cat(1, a1.original_deconv(:, ~isnan(ROIs)), a2.original_deconv(:, ROIs(~isnan(ROIs))));
trials = cat(1, a1.behavior.trials_ts(:), length(a1.behavior.frame_ts) + a2.behavior.trials_ts(:));
trials = discretize(1:length(x), trials);
trials = trials(:);

if ops.rm_mvt
    vel = cat(1, a1.behavior.unit_vel(:), a2.behavior.unit_vel(:));
    thres = noRun(vel);
    thres = (abs(vel) > thres) & ...
        cat(1, ...
        (a1.behavior.trials(1) < a1.behavior.frame_ts(:)) & ...
        (a1.behavior.trials(end) > a1.behavior.frame_ts(:)), ...
        (a2.behavior.trials(1) < a2.behavior.frame_ts(:)) & ...
        (a2.behavior.trials(end) > a2.behavior.frame_ts(:)) ...
        );
    cv_idx = cv_idx(thres);
    x = x(thres);
    n = n(thres, :);
    trials = trials(thres);
end

fs = mean([a1.fs, a2.fs]);
ops.dt = round(fs * ops.dt);
ops.smooth = round(fs * ops.smooth);

md = sam_mape(x, n, trials, 'dt', ops.dt, 'smooth', ops.smooth, 'bins', ops.bins, ...
    'circ', ops.circ, 'cv', true, 'mle', ops.mle, 'penalty', ops.penalty, 'custom', cv_idx);

decoded = md.decoded;
P = md.ll;
pos = md.x;
err = [md.err(:), md.sem(:)];


function ops = parse_inputs(inputs)
ops.dt = 1; % all in seconds
ops.smooth = .2;
ops.bins = 50;
ops.circ = true;
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