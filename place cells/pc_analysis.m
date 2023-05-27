%% Analysis code for place cells
already=0;
vr_length=round(range(unit_pos));
if exist('arms')
    trials=frame_ts(arms);
    trials_ts=arms;
end

%% denoise raw fluorescence signal
% den=medfilt1(tcs.ratio,13,2);
tcs.ratio=double(tcs.ratio);
den=deconv;

sr=1/mean(diff(frame_ts));

%% Plot trials time series
plot(frame_ts,unit_pos); hold on;
plot(frame_ts(trials_ts),unit_pos(trials_ts),'+'); hold off;


%% No-running epochs
[c,centres]=hist((abs(unit_vel)),length(unit_vel)/2);
% [c,centres]=hist(log(abs(unit_vel)),floor(range(abs(unit_vel))/3)); %log-dist
thres=findpeaks(-c);
thres=(centres(thres.loc(1)));
% thres=exp(centres(thres.loc(2))); %log-dist
figure;
hist((abs(unit_vel)),length(unit_vel)/2);
hold on;
plot([thres thres],[0 max(c)],'r');
hold off
figure;
plot(frame_ts,unit_pos); hold on
plot(frame_ts(unit_vel>thres | unit_vel<-thres),unit_pos(unit_vel>thres | unit_vel<-thres),'+');
hold off


%% plot stuff
figure;
ax1=subplot(4,1,[1 2]);
try
    imagesc(Smooth(deconv(:,ordered(end:-1:1))',[0 20]));
    
    %for scatter plot instead
%     plot(peak_rows,peak_cols,'k.','MarkerSize',5);
    
catch
    signal=deconv(:,pc_list);
    imagesc(Smooth(signal(:,ordered(end:-1:1))',[0 5]));
end
set(gca,'xtick',0:2000*sr:length(den));
set(gca,'xticklabel',strsplit(num2str(0:floor(length(den)/sr/2000))));
% xlabel('frames');
ylabel('neuron #');
title(fn);

% subplot(4,1,3);
% plot(1:length(unit_vel),unit_vel);
ax2=subplot(4,1,[3]);
plot(unit_pos);
set(gca,'xtick',0:2000*sr:length(den));
set(gca,'xticklabel',strsplit(num2str(0:floor(length(den)/sr/2000))));
% xlabel('frames');
ylabel('position (cm)');

ax3=subplot(4,1,[4]);
plot(unit_vel);
set(gca,'xtick',0:2000*sr:length(den));
set(gca,'xticklabel',strsplit(num2str(0:floor(length(den)/sr/2000))));
xlabel('time (min)');
ylabel('velocity (cm/s)');

% supertitle({fn});
linkaxes([ax1,ax2,ax3],'x');
colormap hot


%% find per-trial peak responses
ordered_peaks=zeros(trials_ts(1)-1,size(den,2));
for i=1:size(den,2)
    signal(:,i)=smooth(den(:,i),21);
end

for i=1:length(trials_ts)-1
    [~,peaks]=max(signal(trials_ts(i):trials_ts(i+1)-1,:));
    idx=zeros(length(trials_ts(i):trials_ts(i+1))-1,size(den,2));
    for j=1:length(peaks)
        idx(peaks(j),j)=1;
    end
    ordered_peaks=[ordered_peaks; idx];
end
ordered_peaks=[ordered_peaks; zeros(length(trials_ts(end):size(den,1)),size(den,2))];
ordered_peaks=ordered_peaks(:,ordered);

[peak_rows,peak_cols]=find(ordered_peaks);


%% All plots
sd=4;
bins=50;
plotflag=true;

if exist('ordered')
    list=ordered;
elseif exist('pc_list')
    list=pc_list;
else
    list=1:size(deconv,2);
end

if ~already
    [psth,raw_psth,raw_stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,thres,unit_pos,unit_vel,frame_ts,trials);
    already=1;
end

figure;
count=1;
for k=list
    if plotflag
        if count>25
            count=1;
            figure;
        end
        subplot(5,5,count);
        imagesc(psth{k});
        set(gca,'xtick',0:bins/4:bins);
        set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/4:0)));
        title(['n = ' num2str(k)]);
        colormap hot
        ylabel('trials')
        xlabel('distance (cm)')
        colorbar
    end
    count=count+1;
end


%% New stack
if exist('pc_list')
    stack=raw_stack(:,pc_list);
else
    stack=raw_stack;
end
stack=(stack-repmat(min(stack),bins,1));
stack=stack./repmat(max(stack),bins,1);

[~,idx]=max(stack);
[~,ordered]=sort(idx);
stack=stack(:,ordered)';
figure;
imagesc(stack);
set(gca,'xtick',0:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
xlabel('position (cm)');
ylabel('ordered neuron no.');
colormap jet;
c=colorbar; c.Label.String='Norm. Mean dF/F';


%% Q matrix (Battaglia et al.)
% run after making stacks

    qMatrix=corr(stack);

figure;
imagesc(qMatrix);
set(gca,'xtick',0:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
xlabel('position (cm)');
set(gca,'ytick',0:bins/5:bins);
set(gca,'yticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
ylabel('position (cm)');
c=colorbar; c.Label.String='corr. coef.';
colormap jet
axis square

corr_decay=zeros(bins);
corr_decay(1,:)=qMatrix(1,:);
for i=2:bins
    corr_decay(i,:)=[qMatrix(i,1:end-i+1) qMatrix(i,1:i-1)];
end
corr_decay=mean(corr_decay);
corr_decay=[corr_decay corr_decay(end-1:-1:1)];
figure;
plot(-vr_length:(2*vr_length+vr_length/bins)/(bins*2-1):vr_length,corr_decay);
xlabel('distance lag');
ylabel('mean corr');


%% Single-cell sparsity index
sparseness=(sum(stack/bins,2).^2) ./ (sum(stack.^2,2) / bins);
figure
plot(smooth(sparseness,5));
set(gca,'xtick',0:bins/10:bins);
set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/10:0)));
xlabel('neurons');
ylabel('single cell sparsity index');

%% Sparsity
clear sparsity

lamb=raw_stack;
m_lamb=mean(lamb);
Pi=Pi./sum(Pi,2);
Pi=Pi';

sparsity=sum(Pi.*lamb).^2./sum(Pi.*lamb.^2);

try
    sparsity=sparsity(1,pc_list);
catch
end

[c,centres]=hist(sparsity,floor(length(sparsity)/4));
c=cumsum(c./sum(c));

figure;
plot(centres,c);
xlabel('sparsity (%)');
ylabel('cumulative freq.');
axis square


%% Place field width
baseline_thres=range(raw_stack).*.2+min(raw_stack);
width_series=raw_stack>baseline_thres;
try
    width_series=width_series(:,pc_list);
catch
end
for i=1:size(width_series,2)
    temp=width_series(:,i)';
    start=strfind(temp,[0 1]);
    ending=strfind(temp,[1 0]);
    if temp(1)==1
        start=[1 start];
    end
    if temp(end)==1
        ending=[ending length(temp)];
    end
    temp=ending-start;
    width(i)=max(temp);
end
width=width.*vr_length./bins;
[c,centres]=hist(width,floor(length(width)/4));
c=cumsum(c./sum(c));

figure;
plot(centres,c);
xlabel('place field width (cm)');
ylabel('cumulative freq.');
axis square


%% Spatial Tuning Vector (below methods adapted from that Losonzcy paper)
% This is an alternative and perhaps overly complicated method for doing
% the above
% Caveat: The original paper did not work with deconvolved signals, but
% instead marked onsets of transients as epochs with >=2sd from baseline mean
% Here we define onsets as epochs with deconv>0 while disregarding firing
% rates

running_epochs=find(unit_vel<-0.1 | unit_vel>0.1); % velocity thresholds for running epochs
idx=diff(running_epochs);
running_epochs=running_epochs(idx<5); %allow 5 resting epochs
idx=zeros(1,length(den));
idx(running_epochs)=1;
running_epochs=logical(idx);

bins=50;
spatial_vect=zeros(bins,size(den,2));
edges=min(unit_pos(running_epochs)):range(unit_pos(running_epochs))/bins:max(unit_pos(running_epochs));
onsets=den>0;

for i=1:length(edges)-1
    idx=unit_pos>=edges(i) & unit_pos<edges(i+1);
    idx=idx.*running_epochs;
%     spatial_vect(i,:)=sum(onsets.*repmat(idx',1,size(den,2)))./(length(find(idx))/length(find(running_epochs)));
    spatial_vect(i,:)=sum(onsets.*repmat(idx',1,size(den,2)))./length(find(running_epochs));
end

[x,y]=pol2cart(repmat(0:2*pi/bins:2*pi-2*pi/bins,size(den,2),1),spatial_vect');
[spatial_theta,spatial_rho]=cart2pol(sum(x,2),sum(y,2));


%% Tuning vector shuffle test
shuffles=100;

idx=onsets.*repmat(running_epochs',1,size(den,2));
n=arrayfun(@(x) length(find(idx(:,x))),1:size(den,2));

for i=1:shuffles
    shuffle_vect=zeros(bins,size(den,2));
    idx=arrayfun(@(x) ceil(sum(running_epochs).*(rand(1,x))),n,'uniformoutput',false);
    s_onsets=zeros(size(den));
    idx2=find(running_epochs);
    for j=1:length(idx)
        s_onsets(idx2(idx{j}),j)=1;
    end
    for j=1:length(edges)-1
        idx=unit_pos>=edges(j) & unit_pos<edges(j+1);
        idx=idx.*running_epochs;
    %     spatial_vect(i,:)=sum(onsets.*repmat(idx',1,size(den,2)))./(length(find(idx))/length(find(running_epochs)));
        shuffle_vect(j,:)=sum(s_onsets.*repmat(idx',1,size(den,2)))./length(find(running_epochs));
    end
    [x,y]=pol2cart(repmat(0:2*pi/bins:2*pi-2*pi/bins,size(den,2),1),shuffle_vect');
    [shuffle_theta{i},shuffle_rho{i}]=cart2pol(sum(x,2),sum(y,2));
end
for i=1:size(den,2)
    idx=arrayfun(@(x) shuffle_rho{x}(i),1:length(shuffle_rho));
    pval(i)=sum(idx>spatial_rho(i))/shuffles;
end


%% Spatial Info
clear SI

lamb=raw_stack;
m_lamb=mean(lamb);
Pi=Pi./sum(Pi,2);
Pi=Pi';

SI_series=Pi.*lamb./m_lamb.*log2(lamb./m_lamb);
SI=sum(SI_series);

try
    SI=SI(1,pc_list);
catch
end

[c,centres]=hist(SI,floor(length(SI)/4));
c=cumsum(c./sum(c));

figure;
plot(centres,c);
xlabel('spatial info. (bits)');
ylabel('cumulative freq.');
axis square


%% SI shuffle test (circ shift)
useGPU=false;
shuffles=1000;

tic
for i=1:shuffles
    perm=ceil(rand(1)*size(deconv,1));
    shuffled_den=[deconv(perm:end,:);deconv(1:perm-1,:)];
    
    [~,~,lamb,Pi]=getStack(bins,sd,vr_length,shuffled_den,thres,unit_pos,unit_vel,frame_ts,trials);
    m_lamb=mean(lamb);
    Pi=Pi./sum(Pi,2);
    Pi=Pi';
    
    temp=Pi.*lamb./m_lamb.*log2(lamb./m_lamb);
    SI=[SI;sum(temp)];
end
toc

pval=1-sum(SI(1,:)>SI(2:end,:))./shuffles;


%% Bayesian inference decoding
% P(t|d) = P(d|t)P(t)/P(d)
% For now, unique lazy way to get P(d); legit way should be Poisson process
% Taken directly from hd_analysis since the underlying concept is exactly
% the same

theta_bins=25; % number of bins for theta
vel_thres=0.05; % kick out frames with below threshold velocity from prob distributions Pt and Pdt

clear edges

theta_bins=min(unit_pos):range(unit_pos)/theta_bins:max(unit_pos);
theta_idx=arrayfun(@(x) find(unit_pos>=theta_bins(x) & unit_pos<theta_bins(x+1) & unit_vel>vel_thres),1:length(theta_bins)-1,'uniformoutput',false);

Pt=cellfun(@(x) length(x)/length(unit_pos),theta_idx);

Pd=zeros(size(den));
for i=1:size(den,2)
    [pdf,~,idx]=histcounts(den(:,i),1000);
    Pd(:,i)=pdf(idx)./size(den,1);
    for j=1:length(theta_idx)
        [pdf,edges{i,j},idx]=histcounts(den(theta_idx{j},i),50);
        Pdt{i,j}=pdf./length(idx);
    end
end

Ptd=zeros(length(theta_idx),size(den,2));
for i=1:length(unit_pos)
    for j=1:length(theta_idx)
        for k=1:size(den,2)
            idx=find(edges{k,j}>den(i,k),1)-1;
            if isempty(idx)
                Ptd(j,k)=0;
            else
                Ptd(j,k)=Pdt{k,j}(idx)*Pt(j)/Pd(i,k);
            end
        end
    end
    p=prod(Ptd,2);
    [~,idx]=max(p);
    decoded(i)=theta_bins(idx);
%     decoded(i)=sum(p.*idx)./sum(p);
end

plot(unit_pos,'linewidth',1)
hold on
plot(smooth(decoded,21),'linewidth',1)
legend('stimulus','decoded')
set(gca,'xtick',0:60*sr:length(den));
set(gca,'xticklabel',strsplit(num2str(0:floor(length(den)/sr/60))));
xlabel('time (s)')
ylabel('position (cm)')
hold off

%% Template matching
% set variable match (the rest series)
cf_range=1:20;
smoothing=5;
bins=size(raw_stack,1);
trial_len=ceil(mean(diff(trials_ts)));

corr_series=cell(length(cf_range),1);
for cf=cf_range
    if cf>1
        window=arrayfun(@(x) mean(raw_stack((x-1)*cf+1:x*cf,:)),1:floor(size(raw_stack,1)/cf),'uniformoutput',false);
        window=cell2mat(window');
        window=Smooth(window,[smoothing 0]);
        window=reshape(window,1,[]);
    else
        window=Smooth(raw_stack,[smoothing 0]);
        window=reshape(window,1,[]);
    end
    window=flip(window);
    
    for i=1:size(match,1)-floor(size(raw_stack,1)/cf)+1
        target=reshape(Smooth(match(i:floor(size(raw_stack,1)/cf)+i-1,ordered),[smoothing 0]),[],1)';
        corr_series{cf}=[corr_series{cf} corr(window',target')];
    end
end

c_series=zeros(size(deconv,1),cf);
for i=1:length(corr_series)
    c_series(1:length(corr_series{i}),i)=corr_series{i};
end

ax1=subplot(2,1,1);
signal=match-min(match);
signal=signal./max(signal);
imagesc(Smooth(signal(:,ordered)',[0 5]))
colormap jet
title('Rest 2 (post-exposure) Template matching sequence');
title('Rest (post-exposure)');
xlabel('frames');
ylabel('ordered neurons');

cf=17;

ax2=subplot(2,1,2);linkaxes([ax1,ax2],'x');
imagesc(Smooth(c_series,[1 0])')
c=colorbar; c.Label.String='Pearson corr score';
% plot(smooth(corr_series{cf},1))
% title(['Compression factor = ' num2str(trial_len/bins*cf)])
% title(['Compression factor = ' num2str(cf)])
xlabel('frames');
ylabel('cf')
% hold on;
% plot([0 18000],[lol*2+mean(corr_series{17}) lol*2+mean(corr_series{17})],'r--')
% plot([0 18000],[-lol*2-mean(corr_series{17}) -lol*2-mean(corr_series{17})],'r--')
% hold off

figure;
mean_corr=cellfun(@(x) mean(abs(x)),corr_series);
plot(trial_len/bins.*cf_range,mean_corr);
xlabel('compression factor');
ylabel('mean correlation');


%% FR-speed regression
count=1;
for i=1:size(den,2)
    mdl=fitlm(unit_vel,tcs.ratio(:,i));
    if mdl.Coefficients.pValue(2)<0.05
        slope(count)=mdl.Coefficients.Estimate(2);
        rsquared(count)=sign(slope(count))*mdl.Rsquared.Adjusted;
        count=count+1;
    end
end
hist(rsquared,25);
xlabel('Speed vs Firing rate correlation');
ylabel('Proportion cells');


%% Corr with velocity
bins=50;
stack=zeros(bins,size(den,2));
edges=min(unit_vel):range(unit_vel)/bins:max(unit_vel);
for i=1:length(edges)-1
    stack(i,:)=
end

%% Speed vs position
figure;hold on;
% for t=1:size(vel_stack,1)
%     plot(1:bins,vel_stack(t,:),'-','linewidth',1,'color',[0.5 0.5 0.5])
% end

plot(1:bins,nanmean(vel_stack),'k-','linewidth',4)
shadedErrorBar(1:bins,nanmean(vel_stack),nanstd(vel_stack));
set(gca,'xtick',1:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
xlabel('position (cm)');
ylabel('Speed (cm/s)');
hold off

%% Export all figures
name='prelim_data_36';

figHandles = get(groot, 'Children');
for i=1:length(figHandles)
    export_fig(name,'-pdf','-append',figHandles(i));
%     print(figHandles(i),num2str(i),'-dpng');
end





