function corr(obj)
% build correlation matrix

idx=[];
for i=1:length(obj.SCE.on)
    idx=[idx find(obj.SCE.on(i)==obj.ts) : find(obj.SCE.on(i)+obj.SCE.dur(i)==obj.ts)];
end

obj.R=corr(obj.deconv(idx,:));

idx=setxor(idx,1:size(obj.deconv,1));
idx=setxor(idx, find(isnan(sum(obj.deconv,2))));

obj.null_R=corr(obj.deconv(idx,:));