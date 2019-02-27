function filter_bands(obj)

%delta
[a,b,c,d]=butter(2,[1 4]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.delta=filtfilt(sos,g,obj.lfp)./obj.f60_env;

%theta
[a,b,c,d]=butter(2,[5 10]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.theta=filtfilt(sos,g,obj.lfp)./obj.f60_env;

%gamma
[a,b,c,d]=butter(2,[30 140]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.gamma=filtfilt(sos,g,obj.lfp)./obj.f60_env;

%ripples
if obj.fs < 500
    error('LFP has been down sampled too much, the minimum Fs should be 500 Hz');
end
[a,b,c,d]=butter(2,[150 250]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.swr=filtfilt(sos,g,obj.lfp)./obj.f60_env;

% movements are characterized by increase in power in the high frequencies ranging from 500 to 1000 Hz
[a,b,c,d]=butter(2,[500 1000]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp_mvt=filtfilt(sos,g,obj.lfp)./obj.f60_env;

% amplitude envelope for SWR band and ripples detection
env=envelope(obj.swr);
ts=mean(env)+3*std(env);
ts=env>ts;

gaps=.02*obj.fs; %consecutive ripples less than 20 ms apart are concatenated
ts=fill_gaps(ts,gaps);

heads=find(get_head(ts));
tails=get_head(ts(end:-1:1));
tails=find(tails(end:-1:1));

if heads(1)==1
    heads(1)=[];
    tails(1)=[];
end
if tails(end)==length(env)
    heads(end)=[];
    tails(end)=[];
end

obj.swr_env=env;
obj.swr_peaks=zeros(length(heads),1);
obj.swr_on=zeros(length(heads),1);
obj.swr_dur=zeros(length(heads),1);
for i=1:length(heads)
    obj.swr_on(i)=obj.t(heads(i));
    [~,pks]=max(env(heads(i):tails(i)));
    obj.swr_peaks(i)=obj.t(heads(i)+pks-1);
    try
        [~,pks]=findpeaks(env(heads(i):tails(i)));
    catch
        pks=[];
    end
    obj.swr_cyc(i)=length(pks);
    obj.swr_dur(i)=obj.t(tails(i)) - obj.t(heads(i));
end

short=obj.swr_cyc<4; %ripples containing less than 4 cycles are discarded
obj.swr_peaks(short)=[];
obj.swr_on(short)=[];
obj.swr_dur(short)=[];
obj.swr_cyc(short)=[];