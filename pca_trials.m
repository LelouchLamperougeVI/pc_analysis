function trials=pca_trials(pc,detrendFlag)

if nargin>1
    for i=1:size(pc,2)
        pc(:,i)=detrend(pc(:,i));
    end
end

for i=1:size(pc,2)
    sd=std(pc(:,i));
    peaks=findpeaks(pc(:,i));
    peaks=peaks.loc;
    peaks(pc(peaks,i)<2*sd+mean(pc(:,i)))=[];
%     peaks(pc(peaks,i)<2*sd)=[];
    down_peaks=findpeaks(-pc(:,i));
    down_peaks=down_peaks.loc;
    
    trials{i}=[];
    for j=1:length(peaks)
        idx=find(peaks(j)<down_peaks,1);
        try
            trials{i}=[trials{i} [down_peaks(idx-1); down_peaks(idx)]];
        catch
        end
    end
end
