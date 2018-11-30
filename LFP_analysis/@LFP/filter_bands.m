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
[a,b,c,d]=butter(2,[150 200]/obj.fs*2,'bandpass');
[sos,g]=ss2sos(a,b,c,d);
obj.swr=filtfilt(sos,g,obj.lfp)./obj.f60_env;