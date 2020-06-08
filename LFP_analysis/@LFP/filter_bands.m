function filter_bands(obj,filt_60)

if isempty(obj.lfp.lfp)
    return
end

if nargin < 2
    filt_60 = true;
end

if filt_60
    %notch 60
    [a,b,c,d]=butter(2,[58 62]/obj.lfp.fs*2,'stop');
    [sos,g]=ss2sos(a,b,c,d);
    obj.lfp.lfp=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;
    %harmonic
    [a,b,c,d]=butter(2,[178 182]/obj.lfp.fs*2,'stop');
    [sos,g]=ss2sos(a,b,c,d);
    obj.lfp.lfp=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;
end

%delta
[a,b,c,d]=butter(2,[1 4]/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.delta=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%theta
[a,b,c,d]=butter(2,[5 10]/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.theta=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%gamma
[a,b,c,d]=butter(2,[30 140]/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.gamma=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%ripples
if obj.lfp.fs < 500
    error('LFP has been down sampled too much, the minimum Fs should be 500 Hz');
end
[a,b,c,d]=butter(2,[150 250]/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.swr.swr=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

% amplitude envelope for SWR band and ripples detection
env=envelope(obj.lfp.swr.swr);
ts=mean(env)+3*std(env);
ts=env>ts;

gaps=.02*obj.lfp.fs; %consecutive ripples less than 20 ms apart are concatenated
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

obj.lfp.swr.swr_env=env;
obj.lfp.swr.swr_peaks=zeros(length(heads),1);
obj.lfp.swr.swr_on=zeros(length(heads),1);
obj.lfp.swr.swr_dur=zeros(length(heads),1);
for i=1:length(heads)
    obj.lfp.swr.swr_on(i)=obj.lfp.ts(heads(i));
    [~,pks]=max(env(heads(i):tails(i)));
    obj.lfp.swr.swr_peaks(i)=obj.lfp.ts(heads(i)+pks-1);
    try
        [~,pks]=findpeaks(env(heads(i):tails(i)));
    catch
        pks=[];
    end
    obj.lfp.swr.swr_cyc(i)=length(pks);
    obj.lfp.swr.swr_dur(i)=obj.lfp.ts(tails(i)) - obj.lfp.ts(heads(i));
end

short=obj.lfp.swr.swr_cyc < obj.lfp.ops.swr_cyc; %ripples containing less than n cycles are discarded
obj.lfp.swr.swr_peaks(short)=[];
obj.lfp.swr.swr_on(short)=[];
obj.lfp.swr.swr_dur(short)=[];
obj.lfp.swr.swr_cyc(short)=[];

tooClose = [false; diff(obj.lfp.swr.swr_peaks) < obj.lfp.ops.swr_gap];
obj.lfp.swr.swr_peaks(tooClose)=[];
obj.lfp.swr.swr_on(tooClose)=[];
obj.lfp.swr.swr_dur(tooClose)=[];
obj.lfp.swr.swr_cyc(tooClose)=[];