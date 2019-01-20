function sce_spectrum(obj,varargin)
% Extract mean lfp time-freq spectrum of length wdw = [seconds_before seconds_after]
% surrounding each SCE onset, around the frequency limits freq =
% [lower_limit upper_limit]

obj.set_ops(varargin);

freq=obj.ops.freq;
wdw=obj.ops.wdw;
wdw_size=obj.ops.wdw_size;

if isempty(obj.lfp)
    error('no LFP object available');
end

obj.lfp.spectrum(wdw_size,wdw+obj.SCE.on(1));
spec=obj.lfp.spec.spectrum;
for i=2:length(obj.SCE.on)
    obj.lfp.spectrum(diff(wdw),wdw+obj.SCE.on(i));
    spec=spec+obj.lfp.spec.spectrum;
end
spec=spec./length(obj.SCE.on);

% obj.spec=spec;
obj.spec=spec(freq(1)<=obj.lfp.spec.f & freq(2)>=obj.lfp.spec.f, :);