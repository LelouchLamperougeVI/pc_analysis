function [psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack, zscore_stack]=getStack(bins,sd,vr_length,deconv,pos,vel,ft,trials)

sd=sd/vr_length*bins;

zdeconv = zscore(deconv);

raw_psth=zeros(length(trials)-1,bins,size(deconv,2));

trial_bins=discretize(ft,trials);
edges=linspace(min(pos),max(pos),bins+1);
pos_bins=discretize(pos,edges);
if any(isnan(trial_bins)) || any(isnan(pos_bins)); error('yo dun goof''d'); end

Pi = accumarray(pos_bins', 1);
[xx, yy] = ndgrid(pos_bins, 1:size(deconv,2));
raw_stack = accumarray([xx(:) yy(:)], deconv(:)) ./ Pi;
zscore_stack = accumarray([xx(:) yy(:)], zdeconv(:)) ./ Pi;
mu_fr=mean(deconv);

Pt = accumarray([trial_bins' pos_bins'], 1);
for i = 1:size(deconv,2)
    raw_psth(:,:,i) = accumarray([trial_bins' pos_bins'], deconv(:,i)) ./ Pt;
end
vel_stack = accumarray([trial_bins' pos_bins'], vel) ./ Pt;

temp = reshape(permute(raw_psth, [2 1 3]), [bins, (length(trials)-1) * size(deconv,2)]);
temp = fast_smooth(temp, sd);
temp = permute(reshape(temp, [bins, length(trials)-1, size(deconv,2)]), [2 1 3]);
psth = mat2cell(temp, length(trials)-1, bins, ones(1, size(deconv,2)));

stack=fast_smooth(raw_stack,sd);
stack=(stack-min(stack));
stack=stack./max(stack);