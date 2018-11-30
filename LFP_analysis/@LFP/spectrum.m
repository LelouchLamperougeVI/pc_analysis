function spectrum(obj,win)
% generate time-frequency spectrum for visualization
% win: Hamming window length in seconds (default 1 min)

if nargin < 2
    win=60;
end

wdw=win*obj.fs;
nol=floor(.5*wdw);
nfft = max([obj.nfft 2^nextpow2(wdw)]);
f = 0:obj.fs/nfft:obj.fs/2; %estimate to Nyquist

[spec,~,t]=spectrogram(obj.lfp,wdw,nol,nfft,obj.fs);
spec=log(abs(spec));

obj.spec.spectrum = spec;
obj.spec.t = t;
obj.spec.f = f;