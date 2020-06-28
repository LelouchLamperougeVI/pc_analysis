function [psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack]=getStack(vars, deconv, smooth)

if nargin < 3
    smooth = true;
end

sd = vars.sd / vars.vr_length * vars.bins;

% raw_psth=zeros(length(vars.trials)-1,vars.bins,size(deconv,2));

trial_bins=discretize(vars.frame_ts,vars.trials);
edges=linspace(min(vars.unit_pos),max(vars.unit_pos),vars.bins+1);
pos_bins=discretize(vars.unit_pos,edges);
if any(isnan(trial_bins)) || any(isnan(pos_bins)); error('yo dun goof''d'); end

[idx, ~] = find(deconv>0); % new sparse overhaul
temp = deconv(deconv>0);
neur_id = repelem(1:size(deconv, 2), sum(deconv>0, 1));
pos = pos_bins(idx);
trials = trial_bins(idx);
Pt = accumarray([trial_bins' pos_bins'], 1);
raw_psth = accumarray([trials' pos' neur_id'], temp) ./ Pt;
raw_stack = squeeze(sum(raw_psth, 1));

Pi = accumarray(pos_bins', 1);
% [xx, yy] = ndgrid(pos_bins, 1:size(deconv,2));
% temp = deconv(:);
% Pi = accumarray([xx(:) yy(:)], ~isnan(temp));
% temp(isnan(temp)) = 0;
% raw_stack = accumarray([xx(:) yy(:)], temp) ./ Pi;
raw_stack(isnan(raw_stack) | isinf(raw_stack)) = 0;
mu_fr=mean(deconv, 1, 'omitnan');

% Pt = accumarray([trial_bins' pos_bins'], 1);
% for i = 1:size(deconv,2)
%     temp = deconv(:,i);
%     Pt = accumarray([trial_bins' pos_bins'], ~isnan(temp));
%     temp(isnan(temp)) = 0;
%     raw_psth(:,:,i) = accumarray([trial_bins' pos_bins'], temp) ./ Pt;
% end
raw_psth(isnan(raw_psth) | isinf(raw_psth)) = 0;
vel_stack = accumarray([trial_bins' pos_bins'], vars.unit_vel) ./ Pt;

if smooth
    temp = reshape(permute(raw_psth, [2 1 3]), [vars.bins, (length(vars.trials)-1) * size(deconv,2)]);
    temp = fast_smooth(temp, sd);
    temp = permute(reshape(temp, [vars.bins, length(vars.trials)-1, size(deconv,2)]), [2 1 3]);
    psth = mat2cell(temp, length(vars.trials)-1, vars.bins, ones(1, size(deconv,2)));
    
    stack=fast_smooth(raw_stack,sd);
    stack=(stack-min(stack));
    stack=stack./max(stack);
else
    psth = [];
    stack = [];
end