function sce_spectrum(obj,clust,varargin)
% Extract mean lfp time-freq spectrum of length wdw = [seconds_before seconds_after]
% surrounding each SCE onset, around the frequency limits freq =
% [lower_limit upper_limit]

if nargin<2
    clust=0; % compute spectrum for all SCEs irrespective of cluster
end
obj.set_ops(varargin);

freq=obj.ops.freq;
wdw=obj.ops.wdw;
wdw_size=obj.ops.wdw_size;

if isempty(obj.lfp.lfp)
    error('no LFP loaded');
end

if ~clust
    idx=1:length(obj.SCE.on);
else
    idx=find(obj.SCE.clust==clust);
end

obj.spectrum(wdw_size,wdw+obj.SCE.on(idx(1)));
spec_on=obj.spec.spectrum;
for i=2:length(idx)
    obj.spectrum(wdw_size,wdw+obj.SCE.on(idx(i)));
    spec_on=spec_on+obj.spec.spectrum;
end
spec_on=spec_on./length(idx);

obj.spectrum(wdw_size,wdw+obj.SCE.peak(idx(1)));
spec_peaks=obj.spec.spectrum;
for i=2:length(idx)
    obj.spectrum(wdw_size,wdw+obj.SCE.peak(idx(i)));
    spec_peaks=spec_peaks+obj.spec.spectrum;
end
spec_peaks=spec_peaks./length(idx);

rn=10;
% idx=rand(rn,1).*obj.lfp.t(end);
idx=linspace(0,obj.lfp.ts(end)-wdw_size,rn);
obj.spectrum(wdw_size,wdw+idx(1));
null=obj.spec.spectrum;
for i=2:rn
    obj.spectrum(wdw_size,wdw+idx(i));
%     null=null+obj.lfp.spec.spectrum;
    null=[null obj.spec.spectrum];
end
% null=null./rn;

% obj.spec=spec;
obj.spec.spectrum_on=spec_on(freq(1)<=obj.spec.f & freq(2)>=obj.spec.f, :);
obj.spec.spectrum_peak=spec_peaks(freq(1)<=obj.spec.f & freq(2)>=obj.spec.f, :);
obj.spec.f=obj.spec.f(freq(1)<=obj.spec.f & freq(2)>=obj.spec.f);
obj.spec.t=linspace(obj.ops.wdw(1),obj.ops.wdw(2),size(obj.spec.spectrum_on,2));

obj.spec.norm=null(freq(1)<=obj.spec.f & freq(2)>=obj.spec.f, :);