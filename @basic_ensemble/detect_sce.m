function detect_sce(obj,varargin)
%Usage: obj.detect_sce(Name,Value)

obj.set_ops(varargin);

deconv=obj.deconv;
ts=obj.ts;
fs=obj.fs;
sig=obj.ops.sig;
gaps=obj.ops.gaps;

% fs=1/median(diff(ts));
sig=fs*sig;
gaps=round(fs*gaps);


deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
deconv=fast_smooth(deconv,sig);

for ii = 1:length(obj.clust)
    [obj.clust_SCE(ii).SCE, obj.clust_MUA(ii).MUA] = detection(obj.clust{ii});
end
[obj.SCE, obj.MUA] = detection(1:size(deconv,2));


    function [SCE, MUA] = detection(cluster)
        mua=mean(deconv(:,cluster),2,'omitnan');
        thres=mean(mua,'omitnan')+obj.ops.thres*std(mua,'omitnan');
        sce=mua>thres;
        sce=fill_gaps(sce,gaps);
        
        off_thres=mean(mua,'omitnan')+obj.ops.off_thres*std(mua,'omitnan');
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
        
        if length(heads) ~= length(tails)
            alignment=xcorr(heads,tails,5);
            [~,idx] = max(alignment);
            alignment = idx -(5 + 1);
            if length(heads) > length(tails)
                if alignment > 0
                    heads(1:abs(alignment)) = [];
                else
                    heads(end-abs(alignment)+1:end) = [];
                end
            else
                if alignment > 0
                    tails(end-abs(alignment)+1:end) = [];
                else
                    tails(1:abs(alignment)) = [];
                end
            end
        end
        if length(heads) ~= length(tails)
            error('heads and tails are of different lengths; tried aligning but didn''t work');
        end
        
        SCE.dur=ts(tails)-ts(heads);
        SCE.on=ts(heads);
        MUA=mua;
        
        for i=1:length(SCE.on)
            idx=find(ts==SCE.on(i)) : find(ts==SCE.on(i)+SCE.dur(i));
            [~,pk]=max(mua(idx));
            SCE.peak(i)=ts(idx(pk));
        end
        
    end

end
