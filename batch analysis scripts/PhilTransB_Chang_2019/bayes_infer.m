function [decoded, pos, err, err_sem] = bayes_infer(analysis,tau,clust,sig)
% seems like sig=5 gives best accuracy
% bins=80;
bins=50;
sd=4;

deconv=analysis.original_deconv;
behavior=analysis.behavior;

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

vr_length=round(range(unit_pos));
fs=1/median(diff(frame_ts));

tau = tau * fs; % time window length

thres=noRun(unit_vel);

thres=(unit_vel>thres | unit_vel<-thres) & (trials(1) < frame_ts & trials(end) > frame_ts);

artifact = [true diff(unit_vel)~=0];
thres = logical(thres .* artifact);
artifact = [true diff(unit_pos)~=0];
thres = logical(thres .* artifact);

unit_vel=unit_vel(thres);
unit_pos=unit_pos(thres);
frame_ts=frame_ts(thres);
deconv=ca_filt(deconv);
deconv=deconv(thres,:);

deconv=fast_smooth(deconv,sig);
deconv = (deconv - min(deconv)) ./ range(deconv);

trials_idx = knnsearch(frame_ts', trials');

even_trials = trials(2:2:end);
odd_trials = trials(1:2:end);

even_idx=[];
odd_idx=[];
for i=2:2:length(trials)-2
    even_idx = [even_idx trials_idx(i):trials_idx(i+1)];
end
try
    even_idx = [even_idx trials_idx(i+2):trials_idx(i+3)];
    even_trials = [even_trials trials(i+3)];
catch
end
for i=1:2:length(trials)-2
    odd_idx = [odd_idx trials_idx(i):trials_idx(i+1)];
end
try
    odd_idx = [odd_idx trials_idx(i+2):trials_idx(i+3)];
    odd_trials = [odd_trials trials(i+3)];
catch
end
odd_deconv = deconv(odd_idx,:);
odd_unit_pos = unit_pos(odd_idx);
odd_unit_vel = unit_vel(odd_idx);
odd_frame_ts = frame_ts(odd_idx);

% kernel = ones(round(tau),1);
% n = conv2(deconv,kernel,'same');
n = movmean(deconv,round(tau),1);
n = n(even_idx,:);

pos = movmedian(unit_pos,round(tau),1);
pos = pos(even_idx);
pos = discretize(pos, linspace(min(pos), max(pos), bins+1));
% pos = discretize(pos, linspace(-analysis.vr_length, 0, bins+1));

[~,~,stack]=getStack(bins,sd,vr_length,odd_deconv,odd_unit_pos,odd_unit_vel,odd_frame_ts,odd_trials);

decoded = zeros(length(clust)+1, length(pos));
err= zeros(length(clust)+1, bins);
err_sem= zeros(length(clust)+1, bins);
count = 1;
decode(1:length(analysis.psth));
for i = 1:length(clust)
    count=count+1;
    decode(clust{i});
end

    function decode(cluster)
        P = prod(stack(:,cluster).^permute(n(:,cluster),[3 2 1]), 2) .* exp(-tau .* sum(stack(:,cluster),2));
        % P = permute(n,[3 2 1]) .* log(stack);
        % P(isnan(P) | isinf(P)) = 0;
        % P = sum(P, 2) - (tau .* sum(stack,2));
        % P = exp(P);
        P = squeeze(P);
        P = P ./ sum(P,1);
        
        [~,decoded(count,:)] = max(P);
        
        err(count,:) = arrayfun(@(x) mean(min(...
                        [abs(pos(pos == x) - decoded(count, pos == x)) ; ...
                        bins - abs(pos(pos == x) - decoded(count, pos == x))] ...
                        )), 1:bins).* analysis.vr_length/bins;
        err_sem(count,:) = arrayfun(@(x) sem(min(...
                        [abs(pos(pos == x) - decoded(count, pos == x)) ; ...
                        bins - abs(pos(pos == x) - decoded(count, pos == x))] ...
                        ), 2), 1:bins) .* analysis.vr_length/bins;
    end
end