function [psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack]=getStack(bins,sd,vr_length,deconv,pos,vel,ft,trials)

sd=sd/vr_length*bins;

raw_psth=zeros(length(trials)-1,bins,size(deconv,2));

trial_bins=discretize(ft,trials);
edges=linspace(min(pos),max(pos),bins+1);
pos_bins=discretize(pos,edges);
if any(isnan(trial_bins)) || any(isnan(pos_bins)); error('yo dun goof''d'); end

% Pi = accumarray(pos_bins', 1);
[xx, yy] = ndgrid(pos_bins, 1:size(deconv,2));
temp = deconv(:);
Pi = accumarray([xx(:) yy(:)], ~isnan(temp));
temp(isnan(temp)) = 0;
raw_stack = accumarray([xx(:) yy(:)], temp) ./ Pi;
raw_stack(isnan(raw_stack) | isinf(raw_stack)) = 0;
mu_fr=mean(deconv, 1, 'omitnan');

% Pt = accumarray([trial_bins' pos_bins'], 1);
for i = 1:size(deconv,2)
    temp = deconv(:,i);
    Pt = accumarray([trial_bins' pos_bins'], ~isnan(temp));
    temp(isnan(temp)) = 0;
    raw_psth(:,:,i) = accumarray([trial_bins' pos_bins'], temp) ./ Pt;
end
raw_psth(isnan(raw_psth) | isinf(raw_psth)) = 0;
vel_stack = accumarray([trial_bins' pos_bins'], vel) ./ Pt;

temp = reshape(permute(raw_psth, [2 1 3]), [bins, (length(trials)-1) * size(deconv,2)]);
temp = fast_smooth(temp, sd);
temp = permute(reshape(temp, [bins, length(trials)-1, size(deconv,2)]), [2 1 3]);
psth = mat2cell(temp, length(trials)-1, bins, ones(1, size(deconv,2)));

stack=fast_smooth(raw_stack,sd);
stack=(stack-min(stack));
stack=stack./max(stack);