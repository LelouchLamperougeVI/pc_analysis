function [pc_width,pc_loc,M]=ricker_test(signal,psth,frac_trials,sig,plotFlag)
% Unsupervised test for place cells by convolving tuning curve with a
% series of Ricker wavelets of different sigma values.
% Can simultaneously test for significance, place fields width and place
% fields centres.
% This is essentially a continuous wavelet transform technique

if ~exist('plotFlag','var')
    plotFlag=0;
end
if ~exist('sig','var')
    sig=3; %mad threshold; 4 mads correspond rouphly to 6 sd
end

if size(signal,1)~=1
    signal=signal';
end

bins=length(signal);
signal=repmat(signal,1,3);

M=zeros(bins,bins*3);
for n=1:bins
    sd=n/2;
    kernel=ricker_wave(-1.5*bins:1.5*bins,sd);
    M(n,:)=1/sqrt(sd).*conv(signal,kernel,'same');
end
idx=imregionalmax(M);

idx([1 bins],:)=false;

idx=median(median(M))+sig*mad(reshape(M,1,[]),1)<M & idx;
% idx=median(M,2)+sig*mad(M',1)'<M & idx;
idx(:,[1:bins 2*bins+1:3*bins])=[];
M(:,[1:bins 2*bins+1:3*bins])=[];

[pc_width,pc_loc]=find(idx);

if plotFlag
    figure;
    subplot(2,2,1);
    imagesc(M);
    hold on;
    plot(pc_loc,pc_width,'r*');
    ylabel('2*sd');
    subplot(2,2,3);
    plot(signal(bins+1:2*bins));
    xlabel('Location');
    
    subplot(2,2,[2 4]);
    [u,~]=meshgrid(0:2*pi/bins:2*pi,1:size(M,1));
    idx=max(max(M))-min(min(M));
    x=cos(u).*idx+[M M(:,1)].*cos(u);
    y=sin(u).*idx+[M M(:,1)].*sin(u);
    [~,z]=meshgrid(1:size(M,2)+1,1:size(M,1));
    surf(x,y,z,sqrt(x.^2+y.^2));
    axis square
end

pc_width=2.*pc_width;

idx=pc_width>(bins*.8);
pc_loc(idx)=[];
pc_width(idx)=[];

if isempty(pc_width)
    return;
end

for i=1:length(pc_width)
    frac=floor(pc_width(i)/2)-pc_loc(i):floor(pc_width(i)/2)+pc_loc(i);
    frac=mod(frac-1,bins)+1; %circular wrap around indexing thinggy
    frac=psth(:,frac);
    frac=sum(any(frac>0,2))/size(psth,1);
    if frac<frac_trials
        pc_loc(i)=0;
        pc_width(i)=0;
    end
end

idx=pc_loc==0;
pc_loc(idx)=[];
pc_width(idx)=[];