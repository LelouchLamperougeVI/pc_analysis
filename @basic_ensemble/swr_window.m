function [stack,all,out,t]=swr_window(obj)

% wdw=3;
wdw=1;
wdw=floor(obj.lfp.fs_2p * wdw / 2);

t=(-wdw:wdw)/obj.lfp.fs_2p;

peaks=obj.lfp.swr_peaks;
peaks=knnsearch(obj.ts, peaks);
peaks( (peaks+wdw) > length(obj.ts) | (peaks-wdw) < 1 ) = [];

stack=zeros(2*wdw+1, size(obj.deconv,2));
all=[];
idx=[];
count=0;
for i=1:length(peaks)
    temp=obj.deconv(peaks(i)-wdw:peaks(i)+wdw,:);
    if any(isnan(temp(:)))
        continue;
    end
    stack=stack + temp;
    all=[all;temp];
    idx=[idx peaks(i)-wdw:peaks(i)+wdw];
    count=count+1;
end

stack=stack./count;
stack=(stack-mean(obj.deconv,'omitnan'))./std(obj.deconv,'omitnan');

idx=setxor(idx,1:size(obj.deconv,2));
out=obj.deconv(idx,:);

obj.swr_stack=stack;
obj.swr_t=t;