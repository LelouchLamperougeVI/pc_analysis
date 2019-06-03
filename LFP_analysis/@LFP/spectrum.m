function spectrum(obj,win,range)
% generate time-frequency spectrum
% win: Hamming window length in seconds (default 1 min)
% range: 2 elements vector defining the beginning and end of the time
%        window for extraction

if nargin < 2
    win=60;
end
if nargin < 3
    range=[obj.t(1) obj.t(end)];
end

wdw=win*obj.lfp.fs;
nol=floor(.5*wdw);
nfft = max([obj.nfft 2^nextpow2(wdw)]);
f = 0:obj.lfp.fs/nfft:obj.lfp.fs/2; %estimate to Nyquist

lfp=obj.lfp.lfp(obj.lfp.ts>=range(1) & obj.lfp.ts<=range(2));

[spec,~,t]=spectrogram(lfp,wdw,nol,nfft,obj.lfp.fs);
spec=log(abs(spec));

obj.spec.spectrum = spec;
obj.spec.t = t;
obj.spec.f = f;