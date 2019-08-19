function [psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack, zscore_stack]=getStack(bins,sd,vr_length,deconv,pos,vel,ft,trials)

% deconv=deconv(thres,:);

% deconv=log(deconv);
% deconv(isinf(deconv))=nan;
% mean_fr=@(s) exp(mean(s,'omitnan')) .* mean(~isnan(s));

sd=sd/vr_length*bins;
list=1:size(deconv,2);

zdeconv = zscore(deconv);

raw_psth=zeros(length(trials)-1,bins,length(list));
vel_stack=zeros(length(trials)-1,bins);

trial_bins=discretize(ft,trials);
edges=linspace(min(pos),max(pos),bins+1);
pos_bins=discretize(pos,edges);
if any(isnan(trial_bins)) || any(isnan(pos_bins)); error('yo dun goof''d'); end

Pi=zeros(bins,1);
raw_stack=zeros(bins,length(list));
zscore_stack=zeros(bins,length(list));
for i=1:bins
    temp=pos_bins==i;
    raw_stack(i,:)=mean(deconv(temp,:),1);
    zscore_stack(i,:)=mean(zdeconv(temp,:),1);
%     raw_stack(i,:)=mean_fr(deconv(temp,:));
    Pi(i)=sum(temp);
end
raw_stack(isnan(raw_stack))=0;
mu_fr=mean(deconv);
% mu_fr=mean_fr(deconv);

% Pi=zeros(length(trials)-1,bins);
for i=1:length(trials)-1
    temp=trial_bins==i;
    for j=1:bins
        idx=temp & (pos_bins==j);
        raw_psth(i,j,:)=permute(mean(deconv(idx,:),1),[3 1 2]);
%         raw_psth(i,j,:)=mean_fr(deconv(idx,:));
        vel_stack(i,j)=mean(vel(idx),2);
%         Pi(i,j)=sum(idx);
    end
end
raw_psth(isnan(raw_psth) | isinf(raw_psth))=0;

psth=arrayfun(@(x) fast_smooth(raw_psth(:,:,x),sd,2),list,'uniformoutput',false);

% raw_stack=squeeze(mean(raw_psth));
% mu_fr=mean(raw_stack);
% Pi=mean(Pi)';

stack=fast_smooth(raw_stack,sd);
stack=(stack-min(stack));
stack=stack./max(stack);