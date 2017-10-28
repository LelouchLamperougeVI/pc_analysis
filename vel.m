%% infer velocity from position (turns out i didnt even need to do this cuz its already in behavior)
% velocity=diff(behavior.pos_cum)./diff(behavior.ts);
% velocity=smooth(velocity,5);
% plot(behavior.ts(1:end-1),velocity);
% velocity=resample(velocity,mean(diff(behavior.ts)),mean(diff(tcs.tt)));

unit_vel=zeros(1,size(tcs.ratio,1));
unit_pos=zeros(1,size(tcs.ratio,1));
for i=1:size(tcs.ratio,1)
    idx=find(behavior.ts>=tcs.tt(i),1);
%     if idx>length(velocity)
%         idx=idx-1;
%     end
%     unit_vel(i)=velocity(idx);
    unit_pos(i)=behavior.pos_norm(idx);
    unit_vel(i)=behavior.speed(idx);
end

%% denoise raw fluorescence signal
den=medfilt1(tcs.ratio,13,2);

%% plot stuff
subplot(4,1,[1 2]); %raw ratio
% imagesc(tcs.ratio');
imagesc(den');

% h=subplot(4,1,[1 2]); %deconv raster
% % spikeTrain=deconv2spike(deconv,h,1);
% imagesc(deconv');

% subplot(4,1,3);
% plot(1:length(unit_vel),unit_vel);
subplot(4,1,[3 4]);
plot(1:length(unit_pos),unit_pos);
xlim([0 length(unit_pos)-1]);

%% calculate correlation coefficient at 0 lag for all cells (takes a long time)
% R=[];
% progressbar();
% for i=1:size(den,2)
%     parfor j=1:size(den,2)
%         r=xcorr(den(:,[i j]));
%         R=[R r(size(den,2))];
%     end
%     progressbar(i/size(den,2));
% end

high_corr=find(R>=2.5e5);
high_corr=floor(high_corr/size(den,2))+1;
high_corr=unique(high_corr);

%% order sequences using weighted average at positions
vel_thres=0.05;
idx=diff(unit_pos);
idx=find(idx>vel_thres);

avg=sum(den(idx,:));
avg=den(idx,:)./repmat(avg,size(den(idx,:),1),1);
avg=repmat(unit_pos(idx)',1,size(den(idx,:),2)).*avg;
weighted_avg=sum(avg);
[~,ordered]=sort(weighted_avg);
subplot(4,1,[1 2]);
imagesc(den(:,ordered(end:-1:1))');

%% PSTH
global n
global h
global den
global bins
global trials
global frame_ts
global unit_pos
    
n=12;
bins=50;

h=figure;
set(h,'keypressfcn',@psth_plot);

%% Per-trial average ordered sequence
% This is a hastily throw together analysis. Whatever you do please change
% it before publication.
[sorted_pos,idx]=sort(unit_pos);
sorted_den=den(idx,:);
% [~,idx]=max(sorted_den,[],1);
% idx=sorted_pos(idx);
% [~,idx]=sort(idx);
sorted_den=sorted_den(:,ordered(end:-1:1));

bins=100;
ordered_psth=zeros(bins,size(den,2));
bin=range(sorted_pos)/bins;
pre=1;
for i=1:bins
    next=find(sorted_pos>=bin*i+sorted_pos(1),1);
    ordered_psth(i,:)=mean(sorted_den(pre:next,:));
    pre=next;
end

%% For "skip every second" reward trials
global n
global h
global trials_ts
global deconv

sp=mean(unique(diff(tcs.tt))); %sampling period
rp=24; %reward period
i=1;
while(tcs.tt(end)>=(i-1)*rp)
    trials_ts(i)=find(tcs.tt>=(i-1)*rp,1);
    i=i+1;
end

n=1;
h=figure;
set(h,'keypressfcn',@cell_plot);


%% Spatial Tuning Vector (bellow methods adapted from that Losonzcy paper)
% Analysis plan: compute either tuning vector or spatial info score for
% each cell for the running and resting epochs separately. Hypothesis:
% should observe three functional cell types: run-on, run-off and in
% between
bin=100;

sd=std(deconv);
onset=deconv>repmat(2.*sd,length(deconv),1);

%% Spatial Info


%% MRL and Rayleigh test
mrl=zeros(1,size(den,2));
mean_r=mrl;
pval=mrl;
rho=deg2rad(r);
% rho=zeros(length(den),1);
% norm_den=den./repmat(max(den),size(den,1),1); %raw
% norm_den=deconv./repmat(max(deconv),size(deconv,1),1); %deconv
norm_den=den./repmat(sum(den),size(den,1),1);
n=length(den);
for i=1:size(norm_den,2)
    [rx,ry]=pol2cart(rho,norm_den(:,i));
    xy=sum([rx ry]);%./length(rx);
    [mean_r(i),mrl(i)]=cart2pol(xy(1),xy(2));
    z(i)=mrl(i)^2*n;
    pval(i)=exp(sqrt(1+4*n+4*(n^2-mrl(i)^2))-(1+2*n));
end

%% sigh... let's try something else; clearly, something is really wrong
thres=0.01;
denOdeconv=deconv;
denOdeconv(denOdeconv<thres)=0;
n=sum(denOdeconv);
rx=sum(denOdeconv.*repmat(cosd(r),1,size(denOdeconv,2)))./n;
ry=sum(denOdeconv.*repmat(sind(r),1,size(denOdeconv,2)))./n;
mrl=sqrt(rx.^2+ry.^2);
mean_r=cart2pol(rx,ry);
z=mrl.^2.*n;
pval=exp(sqrt(1+4.*n+4*(n.^2-mrl.^2))-(1+2.*n));

%% alright, clearly this shit ain't working. Let's try permutation test
% This test shifts optic angle series from ratio series to compute
% distribution of Kuiper statistic.

delta_bins=80; %divide time series into bins
uni_dist=sorted_rho./2./pi;
[sorted_rho,idx]=sort(rho);
perm_dist=zeros(size(den,1),size(den,2),ceil(size(den,1)/delta_bins));
i=1;
for shift=0:delta_bins:size(den,1)-1
    I=[1:size(den,1)]+shift;
    I(I>size(den,1))=I(I>size(den,1))-size(den,1);
    circ_deconv=deconv(I,:);
%     circ_sorted_deconv=circ_sorted_deconv./repmat(max(circ_sorted_deconv),size(circ_sorted_deconv,1),1);
    perm_dist(:,:,i)=circ_deconv;
    i=i+1;
end

polar_bins=18;
theta_bins=2*pi/polar_bins;
perm.theta=[0:polar_bins-1]*theta_bins;
i=1;
for bin=0:polar_bins-1
    bin_dist=perm_dist(find(rho>=bin*theta_bins & rho<(bin+1)*theta_bins),:,:);
    bin_dist=reshape(bin_dist,size(bin_dist,1)*size(bin_dist,3),size(bin_dist,2));
    bin_mean=mean(perm_dist(find(rho>=bin*theta_bins & rho<(bin+1)*theta_bins),:,:));
    perm.og(i,:)=mean(perm_dist(find(rho>=bin*theta_bins & rho<(bin+1)*theta_bins),:,1));
    perm.mean(i,:)=mean(bin_mean,3);
    perm.median(i,:)=median(bin_dist);
    perm.lowBound(i,:)=prctile(bin_mean,5,3);
    perm.highBound(i,:)=prctile(bin_mean,95,3);
    i=i+1;
end

%% polar plot
n=8;
polarplot([perm.theta perm.theta(1)],perm.lowBound([1:end 1],n),'-.'); hold on;
polarplot([perm.theta perm.theta(1)],perm.highBound([1:end 1],n),'-.');
polarplot([perm.theta perm.theta(1)],perm.mean([1:end 1],n),'-');
polarplot([perm.theta perm.theta(1)],perm.og([1:end 1],n)); hold off;

%%
significance=max(perm.og)'>perm.highBound(repmat(max(perm.og),size(perm.og,1),1)==perm.og);

significance=prctile(V,95);
for i=1:length(significance)
    sig(i)=V(1,i)>significance(i);
end


%% polar plot
clear bin_rho

[sorted_rho,idx]=sort(rho);
circ_sorted_den=den(idx,:);
circ_sorted_deconv=deconv(idx,:);
denOdeconv=circ_sorted_deconv;

bins=50;
bin_theta=0:2*pi/bins:2*pi;
bin_rho_idx=arrayfun(@(x) find(sorted_rho>x,1),bin_theta(1:end-1)');

n=1;
bin_rho(1)=mean(denOdeconv(1:bin_rho_idx(1),n));
for i=1:length(bin_rho_idx)-1
    bin_rho(i+1)=mean(denOdeconv(bin_rho_idx(i):bin_rho_idx(i+1),n));
end
bin_rho(i+2)=mean(denOdeconv(bin_rho_idx(i+1):end,n));
bin_rho=bin_rho./max(bin_rho);
polarplot(bin_theta,bin_rho);hold on;
polarplot([0 mean_r(n)],[0 mrl(n)]);hold off;












