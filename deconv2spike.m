function spikeTrain = deconv2spike(deconv,h,thres)
    if nargin==2
        thres=median(abs(deconv)/0.6745); %autothreshold formula from scholarpedia
    else
        thres=thres.*ones(1,size(deconv,2));
    end
    
    spikeTrain=zeros(size(deconv));
    for i=1:size(deconv,2)
        [~,peaks]=findpeaks(deconv(:,i));
        spikeTrain(peaks,i)=1;
        spikeTrain(deconv(peaks,i)<thres(i),i)=0;
    end
    
    if ishandle(h)
        rasterX=1:size(deconv,1);
        rasterY=spikeTrain.*repmat((1:size(deconv,2)),size(spikeTrain,1),1);
        for i=1:size(rasterY,2)
            hold on;
            plot(h,rasterX,rasterY(:,i),'k.');
        end
    end