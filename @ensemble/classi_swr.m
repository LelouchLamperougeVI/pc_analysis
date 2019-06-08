function classi_swr(obj, wdw, p)
% Classify clusters are either SWR excited (1), inhibited (-1) or not
% modulated by SWRs (0)
%
% Algo:
%   Given two t x N x S matrices: obj.swr_all and obj.swr_null_all
%   Get the mean across all ensemble cluster neurons within a wdw centered
%   around ripple peaks
%   Perform 2 sample one-tailed t-test between sample (mean ensemble 
%   activation measured for S ripples) and null (circshuffled ripples)
%
% Parameters:
%   wdw:    .2 sec  (default)
%   p:      .01     (default)

if nargin<2
    wdw = .2;
end
if nargin<3
    p = .01;
end

if wdw > range(obj.swr_t)
    error('the wdw size exceeds the width of the extraction from ensemble.swr_window');
end

wdw = knnsearch(obj.swr_t', [-wdw; wdw]);

obj.swr_class = zeros(length(obj.clust),1);
for ii = 1:length(obj.clust)
    sample = obj.swr_all(wdw(1):wdw(2), obj.clust{ii}, :);
    denom = size(sample,1) + size(sample,2);
    sample = sum(sample,1);
    sample = sum(sample,2);
    sample = sample ./ denom;
    
    null = obj.swr_null_all(wdw(1):wdw(2), obj.clust{ii}, :);
    denom = size(null,1) + size(null,2);
    null = sum(null,1);
    null = sum(null,2);
    null = null ./ denom;
    
    [~,pval] = ttest2(sample, null, 'tail','right', 'vartype','unequal');
    if pval < p
        obj.swr_class(ii) = 1;
        continue
    end
    [~,pval] = ttest2(sample, null, 'tail','left', 'vartype','unequal');
    if pval < p
        obj.swr_class(ii) = -1;
    end
end