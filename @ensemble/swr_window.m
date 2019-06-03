function [stack,all,out,t]=swr_window(obj)

% wdw=3;
wdw=2;
wdw=floor(obj.twop.fs * wdw / 2);

t=(-wdw:wdw)/obj.twop.fs;

peaks=obj.lfp.swr.swr_peaks;
peaks=knnsearch(obj.twop.ts, peaks);
peaks( (peaks+wdw) > length(obj.twop.ts) | (peaks-wdw) < 1 ) = [];

stack=zeros(2*wdw+1, size(obj.twop.deconv,2));
all=[];
idx=[];
count=0;
for i=1:length(peaks)
    temp=obj.twop.deconv(peaks(i)-wdw:peaks(i)+wdw,:);
    if any(isnan(temp(:)))
        continue;
    end
    stack=stack + temp;
    all=cat(3, all, temp);
    idx=[idx peaks(i)-wdw:peaks(i)+wdw];
    count=count+1;
end

stack=stack./count;
stack=(stack-mean(obj.twop.deconv,'omitnan'))./std(obj.twop.deconv,'omitnan');

idx=setxor(idx,1:size(obj.twop.deconv,2));
out=obj.twop.deconv(idx,:);

obj.swr_stack=stack;
obj.swr_t=t;
obj.swr_all = all;