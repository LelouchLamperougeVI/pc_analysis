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
[a,b,c,d]=butter(2,obj.lfp.ops.freqs('delta')/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.delta=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%theta
[a,b,c,d]=butter(2,obj.lfp.ops.freqs('theta')/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.theta=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%gamma
[a,b,c,d]=butter(2,obj.lfp.ops.freqs('gamma')/obj.lfp.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.lfp.gamma=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;

%ripples
if obj.lfp.fs < 500
    error('LFP has been down sampled too much, the minimum Fs should be 500 Hz');
end
% [a,b,c,d]=butter(2,[150 250]/obj.lfp.fs*2,'bandpass');
% [sos,g]=ss2sos(a,b,c,d);
% obj.lfp.swr.swr=filtfilt(sos,g,obj.lfp.lfp)./obj.lfp.f60_env;
wdw = hamming(obj.lfp.ops.filt_order+1); % creat hamming window
swr_flt = fir1(obj.lfp.ops.filt_order, obj.lfp.ops.freqs('swr')/(obj.lfp.fs/2), 'bandpass', wdw, 'scale');
obj.lfp.swr.swr = filtfilt(swr_flt, 1, obj.lfp.lfp);

try
    % amplitude envelope for SWR band and ripples detection
    wdw = round(obj.lfp.fs * obj.lfp.ops.env_len); % 8ms sliding window
    env=envelope(obj.lfp.swr.swr, wdw, 'rms');
    if ~isempty(obj.behavior)
        mvt_idx = obj.behavior.ts(obj.behavior.speed ~= 0);
        mvt_idx = knnsearch(obj.lfp.ts', mvt_idx');
        mvt_idx = cat(1, mvt_idx, knnsearch(obj.lfp.ts', obj.twop.ts(any(isnan(obj.twop.deconv), 2))));
        mvt = false(length(obj.lfp.ts), 1);
        mvt(mvt_idx) = true;

        kernel = ones(round(obj.lfp.fs), 1); % dilate to 1 sec
        mvt = logical(conv(mvt, kernel, 'same'));

        mvt = fill_gaps(mvt, round(obj.lfp.fs)); % fill 1 sec gaps

        env(mvt) = nan;
    else
        warning('Behavioural data unavailable. Please run LFP.extract_behaviour(). SWR detection will include movement epochs.');
    end
    
    ts = mean(env, 'omitnan') + obj.lfp.ops.swr_thres*std(env, 'omitnan');
    ts = env>ts;
    % identify ripple onsets/offsets
    onset = env > (mean(env, 'omitnan') + obj.lfp.ops.swr_on_thres * obj.lfp.ops.swr_thres * std(env, 'omitnan'));
    heads = find(get_head(onset));
    tails = get_head(onset(end:-1:1));
    tails = find(tails(end:-1:1));
    for ii = 1:length(heads)
        if ~any(ts(heads(ii):tails(ii)))
            onset(heads(ii):tails(ii)) = false;
        end
    end
    ts = onset;

    gaps = obj.lfp.ops.swr_gap * obj.lfp.fs; %consecutive ripples less than 20 ms apart are concatenated
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

%     tooClose = [false; diff(obj.lfp.swr.swr_peaks) < obj.lfp.ops.swr_gap];
%     obj.lfp.swr.swr_peaks(tooClose)=[];
%     obj.lfp.swr.swr_on(tooClose)=[];
%     obj.lfp.swr.swr_dur(tooClose)=[];
%     obj.lfp.swr.swr_cyc(tooClose)=[];
    
    obj.lfp.swr.rate = length(obj.lfp.swr.swr_on)  * obj.lfp.fs / sum(~isnan(env));
    
catch
end