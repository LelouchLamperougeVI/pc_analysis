function reference60(obj,len)
%use 60 Hz as reference to normalize power
%generates signal envelope with Hilbert filter length of wdw (default 1 sec)
if nargin < 2
    len=1;
end
f60=designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5,'SampleRate',obj.lfp.fs);
f60=filtfilt(f60,obj.lfp.lfp);
obj.lfp.f60_env=envelope(f60,floor(obj.lfp.fs*len));