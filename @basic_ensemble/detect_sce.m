function detect_sce(obj,varargin)
%Usage: obj.detect_sce(Name,Value)

obj.set_ops(varargin);

deconv=obj.deconv;
ts=obj.ts;
sig=obj.ops.sig;
thres=obj.ops.thres;
off_thres=obj.ops.off_thres;
gaps=obj.ops.gaps;

fs=1/median(diff(ts));
sig=fs*sig;
gaps=round(fs*gaps);


deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
deconv=fast_smooth(deconv,sig);

mua=mean(deconv,2,'omitnan');
thres=mean(mua,'omitnan')+thres*std(mua,'omitnan');
sce=mua>thres;
sce=fill_gaps(sce,gaps);

off_thres=mean(mua,'omitnan')+off_thres*std(mua,'omitnan');
off=mua<off_thres;
off=off.*~sce;
off=find(off);

heads=get_head(sce);
tails=get_head(sce(end:-1:1));
tails=tails(end:-1:1);
heads=find(heads);tails=find(tails);
if heads(1)==1
    heads(1)=[];
    tails(1)=[];
end
if tails(end)==length(sce)
    heads(end)=[];
    tails(end)=[];
end

idx=heads'>off;
idx=get_head(idx(end:-1:1,:));
idx=idx(end:-1:1,:);
[idx,~]=ind2sub(size(idx),find(idx));
heads=off(idx);

idx=tails'<off;
idx=get_head(idx);
[idx,~]=ind2sub(size(idx),find(idx));
tails=off(idx);

obj.SCE.dur=ts(tails)-ts(heads);
obj.SCE.on=ts(heads);
obj.MUA=mua;

for i=1:length(obj.SCE.on)
    idx=find(ts==obj.SCE.on(i)) : find(ts==obj.SCE.on(i)+obj.SCE.dur(i));
    [~,pk]=max(mua(idx));
    obj.SCE.peak(i)=ts(idx(pk));
end

