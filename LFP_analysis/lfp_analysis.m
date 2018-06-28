function analysis=lfp_analysis(varargin)
% Inputs:
%   lfp,frame_ts,deconv,fs
%   'down_sample':  subsampling rate (default 400 Hz)
%   'win':          stf window size (default 60s)
%   'sig':          number of sd for envelope detection of ripples (default 3.5)
%   'reference60':  use 60Hz to normalize power (default false)
%   'precision':    precision for REM/NREM classification (default 1s)
%   'rta_win':      ripple-triggered averaging window size (default 1s)
%   'r_lag':        minimum lag between SWR events (default 200 ms)
%   'plotFlag':     plot relevant figures (default false)

lfp=varargin{1};
frame_ts=varargin{2};
deconv=varargin{3};
fs=varargin{4};

[down_sample,sig,win,reference60,precision,rta_win,r_lag,plotFlag]=parse_input(varargin(5:end));

ca_fs=1/median(diff(frame_ts));

factor=floor(fs/down_sample);
fs=fs/factor;

lfp=decimate(lfp,factor,'fir');

minute=floor(fs*win);
len=floor(length(lfp)/minute);
data=reshape(lfp(1:len*minute),minute,len);

wdw = 5*down_sample; %5s window
nol = floor(.5*wdw); %half overlapping windows
nfft = 2^16; 
f = 0:fs/nfft:fs/2;
%taking out 60Hz to normalize
fo1=59.5;
fi1=60.5; 
f_60out = f<fo1 | f>fi1;

nblock=size(data,2);
S=zeros(length(f),size(data,2));
for j=1:nblock
    data(:,j)=detrend(data(:,j));
    S(:,j) = pwelch(data(:,j),wdw,nol,nfft,fs);  
    S(:,j)=S(:,j)./sum(S(~f_60out,j)); %use 60hz as reference for power normalization
end
analysis.spec=S;


%%%%%%%%%%
%detrend using 1s window with 500ms steps
lfp=locdetrend(lfp,fs,[1 .5]);

if reference60
    f60=designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5,'SampleRate',fs);
    f60=filtfilt(f60,lfp);
    f60_env=envelope(f60,floor(fs*win));
else
    f60_env=1;
end

rms_t=round(linspace(1,length(lfp),length(lfp)/floor(precision*fs)));

[a,b,c,d]=butter(2,[1.5 4]/fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
delta=filtfilt(sos,g,lfp)./f60_env;
delta_rms=arrayfun(@(x) sqrt(mean(delta(rms_t(x:x+1)).^2)),1:length(rms_t)-1);
% delta_rms=sqrt(mean(reshape(delta,floor(precision*fs),[]).^2));
outliers=isoutlier(delta_rms,'quartiles');
outliers(1)=0;
delta_rms(outliers)=delta_rms(find(outliers)-1);

[a,b,c,d]=butter(2,[5 14]/fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
theta=filtfilt(sos,g,lfp)./f60_env;
theta_rms=arrayfun(@(x) sqrt(mean(theta(rms_t(x:x+1)).^2)),1:length(rms_t)-1);
% theta_rms=sqrt(mean(reshape(theta,floor(precision*fs),[]).^2));
outliers=isoutlier(theta_rms,'quartiles');
outliers(1)=0;
theta_rms(outliers)=theta_rms(find(outliers)-1);

rms_t=(1:length(delta_rms)).*precision;

k=2;
bad_clusters=true;
while(bad_clusters)
    clusters=kmeans([smooth(delta_rms) smooth(theta_rms)],k);
    bad_clusters=abs((sum(clusters==1)-sum(clusters==2))/length(clusters))>.9;
end
% gmm=fitgmdist([delta_rms' theta_rms'],k,'covariancetype','diagonal','options',statset('MaxIter',1000));

analysis.delta=delta;
analysis.delta_rms=delta_rms;
analysis.theta=theta;
analysis.theta_rms=theta_rms;
analysis.rms_t=rms_t;

if plotFlag(1)
    figure;
    ax1=subplot(9, 1, 1:3);
    fo=90;
    fi=150;
    fp = f>fo & f<fi;
    t = (1:nblock).*(win/60);
    pcolor(t,f(fp),10*log10(S(fp,:)))
    shading interp 
    ylabel('frequency (Hz)')

    ax2=subplot(9, 1, 4:6);
    fo=0.5;
    fi=20;
    fp = f>fo & f<fi;
    pcolor(t,f(fp),10*log10(S(fp,:)))
    shading interp 
    ylabel('frequency (Hz)')
    
    ax3=subplot(9, 1, 7);
    t=(1:length(lfp))./fs./60;
    plot(t,theta);
    ylabel('theta band');
    ax4=subplot(9, 1, 8);
    plot(t,delta);
    ylabel('delta band');
    ax5=subplot(9, 1, 9);
    plot(rms_t./60,clusters);
    set(gca,'ytick',1:2,'yticklabel',{'NREM','REM'});
    xlabel('time (min)')
    ylim([0 3]);
    
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
    
    figure;
    hold on
    for i=1:k
        plot(delta_rms(clusters==i),theta_rms(clusters==i),'.');
    end
    xlabel('delta RMS')
    ylabel('theta RMS')
    axis square
end

if mean(delta_rms(clusters==1))>mean(delta_rms(clusters==2))
    nrem=clusters==1;
else
    nrem=clusters==2;
end
nrem(rms_t<frame_ts(1) | rms_t>frame_ts(end))=[];
rem=~nrem;
rms_t(rms_t<frame_ts(1) | rms_t>frame_ts(end))=[];

analysis.rem=rem;
analysis.nrem=nrem;
analysis.rms_t=rms_t;


frame_sample=arrayfun(@(x) find(frame_ts>x,1), rms_t);
ca_fr=zeros(1,length(frame_sample));
for i=1:length(frame_sample)-1
    ca_fr(i)=sum(sum(deconv(frame_sample(i):frame_sample(i+1)-1,:)))./numel(deconv(frame_sample(i):frame_sample(i+1)-1,:));
end
ca_fr(end)=ca_fr(end-1);

gaps=find(get_head(isnan(ca_fr')));
temp=ca_fr(1:gaps(1));
ca_fr(1:gaps(1))=(temp-mean(temp,'omitnan'))./std(temp,'omitnan');
for i=1:length(gaps)-1
    temp=ca_fr(gaps(i):gaps(i+1)-1);
    ca_fr(gaps(i):gaps(i+1)-1)=(temp-mean(temp,'omitnan'))./std(temp,'omitnan');
end
temp=ca_fr(gaps(end):end);
ca_fr(gaps(end):end)=(temp-mean(temp,'omitnan'))./std(temp,'omitnan');
    

rem=ca_fr(rem);
nrem=ca_fr(nrem);

rem_rms=sqrt(mean(reshape(rem(1:50*floor(length(rem)/50)),floor(length(rem)/50),50).^2,'omitnan'));
nrem_rms=sqrt(mean(reshape(nrem(1:50*floor(length(nrem)/50)),floor(length(nrem)/50),50).^2,'omitnan'));

analysis.ca_fr=ca_fr;
analysis.frame_sample=frame_sample;
analysis.rem_rms=rem_rms;
analysis.nrem_rms=nrem_rms;

if plotFlag(2)
    figure;
    ax1=subplot(2,1,1);
    t = (1:nblock).*(win/60);
    fo=0.5;
    fi=20;
    fp = f>fo & f<fi;
    pcolor(t,f(fp),10*log10(S(fp,:)))
    shading interp 
    ylabel('frequency (Hz)')
    ax2=subplot(2, 1, 2);
    plot(rms_t./60,ca_fr)
    ylabel('MUA');
    xlabel('time (min)');
    linkaxes([ax1 ax2],'x');
    
    figure;
    boxplot([rem';nrem'],[zeros(length(rem),1);ones(length(nrem),1)],'labels',{'REM','NREM'})
    ylabel('norm. MUA')
    
    figure;
    boxplot([rem_rms';nrem_rms'],[zeros(length(rem_rms),1);ones(length(nrem_rms),1)],'labels',{'REM','NREM'})
    ylabel('RMS')
end



[a,b,c,d]=butter(2,[90 150]/fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
swr=filtfilt(sos,g,lfp)./f60_env;

swr_env=abs(envelope(swr,ceil(fs/45)));
swr_events=(log(swr_env)-mean(log(swr_env)))>sig*std(log(swr_env));
idx=find(swr_events);
idx([true ; diff(idx)>(r_lag/fs*1000)])=[];
swr_events(idx)=false;

if plotFlag(3)
    figure
    ax1=subplot(4,1,1:2);
    pcolor(frame_ts./60,1:size(deconv,2),fast_smooth(deconv,10)');
    shading interp
    ax2=subplot(4,1,3);
    plot((1:length(swr))./fs./60,swr);
    ax3=subplot(4,1,4);
    plot((1:length(swr))./fs./60,swr_events);
    linkaxes([ax1 ax2 ax3],'x');
end

swr_events=find(swr_events)./fs;

ripple_frames=knnsearch(frame_ts,swr_events);

rta_win=floor(ca_fs*rta_win/2);
rta=zeros(rta_win*2+1,length(ripple_frames));
for i=1:size(rta,2)
    temp=deconv(ripple_frames(i)-rta_win:ripple_frames(i)+rta_win,:);
    rta(:,i)=sum(temp,2)./ca_fs;
end

outliers=isoutlier(sum(rta));
rta(:,outliers)=[];

analysis.swr=swr;
analysis.swr_events=swr_events;
analysis.ripple_frames=ripple_frames;
analysis.rta=rta;

if plotFlag(4)
    figure;
    hold on
    plot(lfp);
    plot(swr-max(swr)-range(lfp));
    plot(range(swr).*swr_events-max(swr)-range(lfp))
    plot(theta-max(theta)-range(swr)-range(lfp))
    plot(delta-max(delta)-range(theta)-range(swr)-range(lfp))
    
    figure
    t=(-rta_win:rta_win)./ca_fs.*1000;
    n=1:size(rta,2);
    pcolor(t,n,fast_smooth(zscore(rta),1)');
    shading interp
    xlabel('time from ripple (ms)')
    ylabel('ripple events')
    
    figure
    shadedErrorBar(t,rta',{@mean,@std})
    xlabel('time from ripple (ms)')
    ylabel('MUA (spk/s)')
end


function [down_sample,sig,win,reference60,precision,rta_win,r_lag,plotFlag]=parse_input(inputs)
plotFlag=false;
sig=3.5;
down_sample=400;
win=60;
reference60=false;
precision=1;
rta_win=1;
r_lag=200;

idx=1;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 'down_sample'
            idx=idx+1;
            down_sample=inputs{idx};
        case 'win'
            idx=idx+1;
            win=inputs{idx};
        case 'reference60'
            idx=idx+1;
            reference60=inputs{idx};
        case 'precision'
            idx=idx+1;
            precision=inputs{idx};
        case 'plotflag'
            idx=idx+1;
            plotFlag=inputs{idx};
        case 'rta_win'
            idx=idx+1;
            rta_win=inputs{idx};
        case 'r_lag'
            idx=idx+1;
            r_lag=inputs{idx};
        otherwise
    end
    idx=idx+1;
end

