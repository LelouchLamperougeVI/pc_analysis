function [stack,all,t]=swr_window(obj)

% wdw=3;
shuffles = 100;
wdw=2;
wdw=floor(obj.twop.fs * wdw / 2);

ts = obj.twop.ts;
deconv = obj.twop.deconv;

t=(-wdw:wdw)/obj.twop.fs;

peaks=obj.lfp.swr.swr_peaks;
peaks=knnsearch(ts, peaks);

[stack, all] = get_stack(peaks, wdw, ts, deconv);

index = ( 1:size(deconv,1) )';
null_all = cell(1,1, shuffles);
xx = size(obj.twop.deconv,1);
parfor ii = 1:shuffles
    idx = circshift(index, randi(xx, 1));
    pks = idx(peaks);
    [~, null_all{ii}] = get_stack(pks, wdw, ts, deconv);
end
null_all = cell2mat(null_all);

obj.ensembles.swr.stack=stack;
obj.ensembles.swr.t=t;
obj.ensembles.swr.all = all;
obj.ensembles.swr.null_all = null_all;



function [stack, all] = get_stack(peaks, wdw, ts, deconv)
peaks( (peaks+wdw) > length(ts) | (peaks-wdw) < 1 ) = [];
stack=zeros(2*wdw+1, size(deconv,2));
all=[];
count=0;
for i=1:length(peaks)
    temp=deconv(peaks(i)-wdw:peaks(i)+wdw,:);
    if any(isnan(temp(:)))
        continue;
    end
    stack=stack + temp;
    all=cat(3, all, temp);
    count=count+1;
end
stack=stack./count;
stack=(stack-mean(deconv,'omitnan'))./std(deconv,'omitnan');

all = (all-mean(deconv,'omitnan'))./std(deconv,'omitnan');

