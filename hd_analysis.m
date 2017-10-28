%% Analysis code for hd cells
den=deconv;

%% Plot trials time series
figure
hold on;
idx=arrayfun(@(x) find(frame_ts>=x,1),trials);
% plot(frame_ts(idx),deg2rad(unit_rot(idx)),'bv','linewidth',3);
% plot(frame_ts,deg2rad(unit_rot),'linewidth',3);
plot(frame_ts,unit_rot,'linewidth',3);
plot(frame_ts(idx),(max(unit_rot)+20).*ones(1,length(idx)),'v');
plot(frame_ts(trials_ts),unit_rot(trials_ts)+5,'v'); 
hold off;
xlabel('time (sec)');
ylabel('rotation (rad)');
legend('reward')

%% Plot Deconv series
signal=Smooth(deconv(:,ordered),[5 0]);
signal=signal-min(signal);
signal=signal./max(signal);
imagesc(signal');
colormap hot
xlabel('frames')
ylabel('ordered neuron no.')

%% All plots
n=97;
skip=1;
bins=90;
quadrants=3;

theta=0:-2*pi/bins:-(2*pi-2*pi/bins);
edges=min(unit_rot):range(unit_rot)/bins:max(unit_rot);
psth=zeros(bins,length(trials_ts)-1,size(deconv,2));
for i=1:length(trials_ts)-1
    for j=1:bins
        idx=find(unit_rot>=edges(j) & unit_rot<edges(j+1));
        idx(idx<trials_ts(i) | idx>=trials_ts(i+1))=[];
        psth(j,i,:)=reshape(mean(deconv(idx,:)),1,1,size(deconv,2));
    end
end
per_trial=psth;
psth=mean(psth,2);
psth=reshape(psth,size(psth,1),size(psth,3));
raw_stack=psth;
psth=Smooth(psth,[2 0]);
psth=(psth-min(psth))./max(psth);
[~,ordered]=max(psth);
[~,ordered]=sort(ordered);

if isempty(n)
    figure;
    imagesc(psth(:,ordered)');
    set(gca,'xtick',0:bins/10:bins);
    set(gca,'xticklabel',strsplit(num2str(0:-2*pi/10:-(2*pi-2*pi/10))));
    xlabel('rotation (rad)');
    ylabel('ordered neurons');
    colormap hot
end

figure;
count=1;
if ~isempty(n)
    subplot(2,1,1);
    imagesc(Smooth(per_trial(:,:,n),[2 0])');
    set(gca,'xtick',0:bins/10:bins);
    set(gca,'xticklabel',strsplit(num2str(0:-2*pi/10:-(2*pi-2*pi/10))));
    xlabel('rotation (rad)');
    ylabel('trials');
    title(['n = ' num2str(n)]);
    colormap hot
    subplot(2,1,2);
    polarplot(theta,psth(:,n));
else
    for i=ordered
        if count>quadrants^2*2
            count=1;
            figure;
        end
        subplot(quadrants,2*quadrants,count);
        imagesc(Smooth(per_trial(:,:,i),[1 0])');
        set(gca,'xtick',0:bins/10:bins);
        set(gca,'xticklabel',strsplit(num2str(0:-2*pi/10:-(2*pi-2*pi/10))));
        xlabel('rotation (rad)');
        ylabel('trials');
        title(['n = ' num2str(i)]);
        colormap hot
        count=count+1;
        subplot(quadrants,2*quadrants,count);
        polarplot(theta,psth(:,i));
        count=count+1;
    end
end


%% All plots with skip
n=[];
bins=80;
quadrants=3;

theta=0:-2*pi/bins:-(2*pi-2*pi/bins);
edges=min(unit_rot):range(unit_rot)/bins:max(unit_rot);
psth=zeros(bins,length(trials_ts)-1,size(deconv,2));
for i=1:length(trials_ts)-1
    for j=1:bins
        idx=find(unit_rot>=edges(j) & unit_rot<edges(j+1));
        idx(idx<trials_ts(i) | idx>=trials_ts(i+1))=[];
        psth(j,i,:)=reshape(mean(deconv(idx,:)),1,1,size(deconv,2));
    end
end

idx=arrayfun(@(x) find(frame_ts>=x,1),trials);
idx=arrayfun(@(x) find(trials_ts>=x,1),idx,'uniformoutput',false);
idx=cell2mat(idx);
idx1=zeros(1,length(trials_ts)-1);
idx1(idx)=1;
idx1=logical(idx1);

per_trial=psth;
psth1=mean(psth(:,idx1,:),2);
psth2=mean(psth(:,~idx1,:),2);
psth1=reshape(psth1,size(psth1,1),size(psth1,3));
psth2=reshape(psth2,size(psth2,1),size(psth2,3));
% raw_stack=psth;
psth1=Smooth(psth1,[2 0]);
psth2=Smooth(psth2,[2 0]);
psth1=(psth1-min(psth1))./max(psth1);
psth2=(psth2-min(psth2))./max(psth2);
[~,ordered]=max(psth);
[~,ordered]=sort(ordered);

figure;
count=1;
if ~isempty(n)
    subplot(2,2,1);
    imagesc(Smooth(per_trial(:,idx1,n),[2 0])');
    set(gca,'xtick',0:bins/6:bins);
    set(gca,'xticklabel',strsplit(num2str(0:-2*pi/6:-(2*pi-2*pi/6))));
    xlabel('rotation (rad)');
    ylabel('trials');
    title('no reward');
    colormap hot
    subplot(2,2,2);
    imagesc(Smooth(per_trial(:,~idx1,n),[2 0])');
    set(gca,'xtick',0:bins/6:bins);
    set(gca,'xticklabel',strsplit(num2str(0:-2*pi/6:-(2*pi-2*pi/6))));
    xlabel('rotation (rad)');
    ylabel('trials');
    title('reward');
    colormap hot
    subplot(2,2,3);
    polarplot(theta,psth1(:,n));
    subplot(2,2,4);
    polarplot(theta,psth2(:,n));
else
    for i=ordered
        if count>quadrants^2*2
            count=1;
            figure;
        end
        subplot(quadrants,2*quadrants,count);
        imagesc(Smooth(per_trial(:,:,i),[1 0])');
        set(gca,'xtick',0:bins/10:bins);
        set(gca,'xticklabel',strsplit(num2str(0:-2*pi/10:-(2*pi-2*pi/10))));
        xlabel('rotation (rad)');
        ylabel('trials');
        title(['n = ' num2str(i)]);
        colormap hot
        count=count+1;
        subplot(quadrants,2*quadrants,count);
        polarplot(theta,psth(:,i));
        count=count+1;
    end
end


%% Plot per-cell trials
global n
global h
global frame_ts
global trials_ts
global trials
global deconv
global den
global mean_vect_rho
global mean_vect_theta
global unit_rot
global pval
global fn
% 
% [~,trials_ts]=findpeaks(unit_rot);
% trials_ts(trials_ts<330)=[];
% % trials_ts=trials_ts.loc;
% trials_ts(1)=[];

n=0;
h=figure;
set(h,'keypressfcn',@hd_plot);
% set(h,'keypressfcn',@hd_plot_noskip);
% set(h,'keypressfcn',@tc_plot);


%% Mean Vector Permutation test
% This test shifts optic angle series from ratio series to compute null
% distribution of mean vector length.
% Like the Rayleigh test, this test does not take into account
% bidirectional tuning.

shuffles=100; %divide time series into number of bins
% perm_dist=zeros(size(den,1),size(den,2),shift_delta);
theta=deg2rad(repmat(unit_rot',1,size(den,2)));
mean_vect_x=zeros(shuffles,size(den,2)); mean_vect_y=mean_vect_x;
for i=1:shuffles-1
    [x,y]=pol2cart(theta,den);
    mean_vect_x(i,:)=sum(x)./size(den,1);
    mean_vect_y(i,:)=sum(y)./size(den,1);
    
    shift=i*length(unit_rot)/shuffles;
    unit_rot2=[unit_rot(shift:end) unit_rot(1:shift-1)];
    theta=deg2rad(repmat(unit_rot2',1,size(den,2)));
end

[mean_vect_theta,mean_vect_rho]=cart2pol(mean_vect_x,mean_vect_y);

% Compute percentile
nless=arrayfun(@(x) mean_vect_rho(:,x)<mean_vect_rho(1,x), 1:size(mean_vect_rho,2),'UniformOutput',false);
nless=cellfun(@(x) sum(x),nless);
nequal=arrayfun(@(x) mean_vect_rho(:,x)==mean_vect_rho(1,x), 1:size(mean_vect_rho,2),'UniformOutput',false);
nequal=cellfun(@(x) sum(x),nequal);
pval=1-(nless + 0.5.*nequal)./size(mean_vect_rho,1);


%% Mean Vector Shuffle Test
% Randomize optic flow patterns
% Not a very legit test, but doing it anyway cuz I'm too dumb to realize I
% shouldn't have used constant velocity

shuffles=200; %number of randomizations
theta=deg2rad(repmat(unit_rot',1,size(den,2)));
mean_vect_x=zeros(shuffles,size(den,2)); mean_vect_y=mean_vect_x;
for i=1:shuffles
    [x,y]=pol2cart(theta,den);
    mean_vect_x(i,:)=sum(x)./size(den,1);
    mean_vect_y(i,:)=sum(y)./size(den,1);
    
    unit_rot2=genRandomRot(frame_ts);
    theta=deg2rad(repmat(unit_rot2',1,size(den,2)));
end

[mean_vect_theta,mean_vect_rho]=cart2pol(mean_vect_x,mean_vect_y);

% Compute percentile
nless=arrayfun(@(x) mean_vect_rho(:,x)<mean_vect_rho(1,x), 1:size(mean_vect_rho,2),'UniformOutput',false);
nless=cellfun(@(x) sum(x),nless);
nequal=arrayfun(@(x) mean_vect_rho(:,x)==mean_vect_rho(1,x), 1:size(mean_vect_rho,2),'UniformOutput',false);
nequal=cellfun(@(x) sum(x),nequal);
pval=1-(nless + 0.5.*nequal)./size(mean_vect_rho,1);


%% Decoding
% Bayesian inference decoding
% P(t|d) = P(d|t)P(t)/P(d)
% For now, unique lazy way to get P(d); legit way should be Poisson process
clear edges;

theta_bins=30; % number of bins for theta

theta_bins=0:360/theta_bins:360-360/theta_bins;
theta_idx=arrayfun(@(x) find(unit_rot>=theta_bins(x) & unit_rot<theta_bins(x+1)),1:length(theta_bins)-1,'uniformoutput',false);
Pt=cellfun(@(x) length(x)/length(unit_rot),theta_idx);

Pd=zeros(size(den));
for i=1:size(den,2)
    [pdf,~,idx]=histcounts(den(:,i),1000);
    Pd(:,i)=pdf(idx)./size(den,1);
    for j=1:length(theta_idx)
        [pdf,edges{i,j},idx]=histcounts(den(theta_idx{j},i),10);
        Pdt{i,j}=pdf./length(idx);
    end
end

Ptd=zeros(length(theta_idx),size(den,2));
for i=1:length(unit_rot)
    for j=1:length(theta_idx)
        for k=1:size(den,2)
            idx=find(edges{k,j}>den(i,k),1);-1;
            if idx==11
                idx=idx-1;
            end
            if isempty(idx)
                Ptd(j,k)=0;
            else
                Ptd(j,k)=Pdt{k,j}(idx)*Pt(j)/Pd(i,k);
            end
        end
    end
    Ptd(:,~all(Ptd))=[]; % cells with any p=0 are excluded; should find a better way to model pdf
    p=prod(Ptd.*10,2);
    [~,idx]=max(p);
    decoded(i)=theta_bins(idx);
%     decoded(i)=sum(p.*idx)./sum(p);
end

plot(unit_rot,'linewidth',1)
hold on
plot(smooth(decoded,11),'linewidth',1)
legend('stimulus','decoded')
xlabel('frames')
ylabel('degrees')
hold off


%% polar plot
clear bin_rho

[sorted_theta,idx]=sort(unit_rot);
circ_sorted_den=den(idx,:);
circ_sorted_deconv=deconv(idx,:);
denOdeconv=circ_sorted_deconv;

bins=30;
bin_theta=0:2*pi/bins:2*pi;
bin_theta_idx=arrayfun(@(x) find(sorted_theta>x,1),bin_theta(1:end-1)');

n=62;
bin_rho(1)=mean(denOdeconv(1:bin_theta_idx(1),n));
for i=1:length(bin_theta_idx)-1
    bin_rho(i+1)=mean(denOdeconv(bin_theta_idx(i):bin_theta_idx(i+1),n));
end
bin_rho(i+2)=mean(denOdeconv(bin_theta_idx(i+1):end,n));
bin_rho=bin_rho./max(bin_rho);
polarplot(bin_theta,bin_rho);hold on;
polarplot([0 mean_vect_theta(1,n)],[0 mean_vect_rho(1,n)]);hold off;