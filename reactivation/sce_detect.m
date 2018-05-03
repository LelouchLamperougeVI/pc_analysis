function [sce,sce_length]=sce_detect(varargin)
% Detect the locations of SCEs using a CWT method
% Description: Compute CWT of deconv MUA activity with a maxican hat mother
% wavelet. Identify local maxima and overlap with significance threshold
% mask (same as in ricker_test). Noise is estimated from a sliding window
% instead of the whole thing.
% Caveat: you should preprocess the deconv for better results (smooth,
% zscore, etc.)
% Inputs:
%   deconv: deconvolved signal
%   'sig': significance threshold in MADs (4 MAD ~= 6 SDs)
%   't_limit': estimated maximum length of SCE event
%   'favor': favor either 'long' or 'short' SCEs (default long)
%   'win': size of moving noise window (default inf)
%   'plot': plot relevant figures (default false)
%   'mem_lim': memory limit (default 6 gb - for my shitty ass computer)
% Output:
%   sce: logical index of SCEs
%   sce_length: duration of each SCE

[sig,t_limit,favor,win,plotFlag,mem_lim]=parse_input(varargin);
deconv=varargin{1};

signal=sum(deconv,2);
wt=cwt_mexh(signal,t_limit)';

maxi=imregionalmax(wt);
if isinf(win)
    thres_mask=sig.*mad(wt,1)<wt;
else
    thres_mask=false(size(wt));
    if win>size(wt,1)
        error('your window is bigger than your recording dumbass -.-''');
    end
    win=floor(win)/2;
    
    thres_mask(1:win+1,:)=wt(1:win+1,:)>sig.*mad(wt(1:win*2+1,:),1);
    
    chunk=floor(mem_lim*1000000000/size(wt,2)/(win*2+1)/8);
    idx=[win+2:chunk:size(wt,1)-win-1 size(wt,1)-win-1];
    for i=1:length(idx)-1
        buff=zeros(win*2+1,(idx(i+1)-idx(i))*size(wt,2));
        for j=1:(idx(i+1)-idx(i))
            buff(:,(j-1)*size(wt,2)+1:j*size(wt,2))=wt(idx(i)+j-win:idx(i)+j+win,:);
        end
        buff=sig.*mad(buff,1);
        buff=reshape(buff,size(wt,2),idx(i+1)-idx(i))';
        thres_mask(idx(i):idx(i+1)-1,:)=wt(idx(i):idx(i+1)-1,:)>buff;
    end
    
%     thres_mask(end-win:end,:)=wt(end-win:end,:)>sig.*mad(wt(end-win*2:end,:),1);
end

sce=maxi&thres_mask;
sce(:,[1 end])=false;
sce=sum(sce,2);

if plotFlag
    figure;
    imagesc(deconv');
    set(gca,'YDir','normal');
    hold on; 
    plot((signal-min(signal))./range(signal).*size(deconv,2)./4,'r'); 
    plot(sce.*size(deconv,2)./4+size(deconv,2)/2,'r');
end


function mads=estim_mad(A,init,margin)
% A more efficient way to estimate local median absolute deviation
idx=1:floor(size(A,1)/init):size(A,1);
mads=zeros(1,length(idx));
for i=1:length(idx)-1
    mads(i)=mad(A(idx(i):idx(i+1)-1),1);
end
mads(end)=mad(A(idx(end):end),1);





function [sig,t_limit,favor,win,plotFlag,mem_lim]=parse_input(inputs)
sig=3;
t_limit=floor(size(inputs{1},1));
favor='long';
win=inf;
plotFlag=false;
mem_lim=6; %in gigabytes

idx=2;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 't_limit'
            idx=idx+1;
            t_limit=inputs{idx};
        case 'favor'
            idx=idx+1;
            favor=inputs{idx};
        case 'win'
            idx=idx+1;
            win=inputs{idx};
        case 'plot'
            idx=idx+1;
            plotFlag=inputs{idx};
        case 'mem_lim'
            idx=idx+1;
            mem_lim=inputs{idx};
        otherwise
    end
    idx=idx+1;
end
