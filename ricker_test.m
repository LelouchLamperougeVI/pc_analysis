function [pc_width,pc_loc,pval]=ricker_test(signal,plotFlag)
% Unsupervised test for place cells by convolving tuning curve with a 
% series of Ricker wavelets of different sigma values.
% Can simultaneously test for significance, place fields width and place
% fields centres.

if nargin<2
    plotFlag=0;
end

if size(signal,1)~=1
    signal=signal';
end

bins=length(signal);

signal=repmat(signal,1,3);
% signal=zscore(signal);
M=zeros(bins);

for n=1:bins
    sd=n/2;
%     alpha=(n-1)/2/sd;
%     kernel=gausswin(n,alpha);
    kernel=ricker_wave(-bins:bins,sd);
%     kernel=zscore(kernel);
    
    temp=conv(signal,kernel,'same');
    M(n,:)=temp(bins+1:2*bins);
end

idx=imregionalmax(M);
[pc_width,pc_loc]=find(idx);
maxi=M(idx);
dist=reshape(M,1,[]);
pval=zeros(1,length(maxi));
for i=1:length(maxi)
    pval(i)=1-sum(maxi(i)>dist)/length(dist);
end

if plotFlag
    figure;
    subplot(2,1,1);
    imagesc(M);
    hold on;
    plot(pc_loc(pval<0.05),pc_width(pval<0.05),'r*');
    subplot(2,1,2);
    plot(signal(bins+1:2*bins));
end