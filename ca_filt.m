function deconv=ca_filt(deconv,plotFlag)
% Filter out left component from calcium signal.
% At this point, I'm not sure what that component is doing there although
% it appears to contain some spatial information.
% Only useful for certain analyses.

if nargin<2
    plotFlag=0;
end

for i=1:size(deconv,2)
    signal=deconv(:,i);
    signal=log(signal);
    idx=isinf(signal);
    signal(idx)=[];
    clusters=kmeans(signal,2);
    clusters=logical(clusters-1);
    if mean(signal(clusters))<mean(signal(~clusters))
        clusters=~clusters;
    end
    
    if any(i==1:3) && plotFlag % plot 3 examples
        figure;
        signal=deconv(:,i);
        signal(signal<=0)=[];
        histlog(signal(clusters),50);
        hold on;
        histlog(signal(~clusters),50);
    end
    
    idx=find(~idx);
    idx(clusters)=[];
    deconv(idx,i)=0;
end