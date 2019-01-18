function [pc_width,pc_loc,M,reject]=ricker_test(signal,psth,frac_trials,sig,width,io_ratio,consecutive,plotFlag)
% Unsupervised test for place cells by convolving tuning curve with a
% series of Ricker wavelets of different sigma values.
% Can simultaneously test for significance, place fields width and place
% fields centres.
% This is essentially a continuous wavelet transform technique
% Inputs:
%   signal:         mean fr vs pos
%   psth:           trials x pos matrix
%   frac_trials:    minimum fraction of trials that need to contain spike (default 1/3)
%   sig:            mad threshold - 4 mads roughly equates to 6 SDs (default 3)
%   width:          minimum and maximum width of place fields (default [.05 .8] of entire track)
%   io_ratio:       ratio between fr inside vs outside of place field (default 2.5)
%   consecutive:    whether the fraction of trials need to be in consecutive trials (default true)
%   plotFlag:       plot pretty figures

if ~exist('plotFlag','var')
    plotFlag=0;
end
if ~exist('sig','var')
    sig=3; %mad threshold; 4 mads correspond rouphly to 6 sd
end
if ~exist('width','var')
    width=[.05 .8];
end
if ~exist('io_ratio','var')
    io_ratio=2.5;
end
if ~exist('frac_trials','var')
    frac_trials=1/3;
end
if ~exist('consecutive','var')
    consecutive=true;
end

if size(signal,1)~=1
    signal=signal';
end

reject='';

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
idx([1:ceil(bins*width(1)/2) floor(bins*width(2)/2):end],:)=false;

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

if isempty(pc_width)
    reject='empty loc_max';
    return;
end

%remove overlapping place fields
[pc_width,idx]=sort(pc_width,'descend');
pc_loc=pc_loc(idx);
for i=1:length(pc_width)
    frac=pc_loc(i)-floor(pc_width(i)/2):floor(pc_width(i)/2)+pc_loc(i);
    frac=mod(frac-1,bins)+1; %circular wrap around indexing thinggy
    idx=ismember(pc_loc,frac);
    idx(i)=false;
    pc_loc(idx)=0;
    pc_width(idx)=0;
end
idx=pc_loc==0;
pc_loc(idx)=[];
pc_width(idx)=[];

% psth segments outside of place fields
out_field=[];
get_outfield;
if length(out_field)<bins*width(1)
    pc_loc=[];
    pc_width=[];
    reject='no out';
    return;
end

% count=length(pc_width);
for i=length(pc_width):-1:1
% while(count)
    get_outfield;
    out_field=psth(:,out_field);
    
    frac=pc_loc(i)-floor(pc_width(i)/2):floor(pc_width(i)/2)+pc_loc(i);
    frac=mod(frac-1,bins)+1; %circular wrap around indexing thinggy
    frac=psth(:,frac);
    %     if (mean_fr(frac) / mean_fr(out_field)) < io_ratio
    if max(mean(frac)) / mean(mean(out_field)) < io_ratio %yes, this is not a mistake, I'm cheating here...
        %     if (sum(frac(frac>prc))/numel(frac)) / (sum(out_field(out_field>prc))/numel(out_field)) < io_ratio
        pc_loc(i)=nan;
        pc_width(i)=nan;
        reject=[reject 'io_rate;'];
        continue;
    end
    
    %     frac=sum(any(frac==max(psth,[],2),2))/size(frac,1);
    frac=any(frac==max([frac out_field],[],2),2) - any(out_field==max([frac out_field],[],2),2);
    frac=frac==1;
    if consecutive
        frac=fill_gaps(frac,max([floor(size(psth,1)*.05) 1])); %allow 5% trials gaps, let's not be so harsh
        frac=find(flipud(get_head(flipud(frac)))) - find(get_head(frac)) + 1;
        frac=max(frac);
    end
    frac=sum(frac)/size(psth,1);
    if frac<frac_trials
        pc_loc(i)=nan;
        pc_width(i)=nan;
        reject=[reject 'frac_trials;'];
    end
end

    function get_outfield()
        out_field=[];
        for k=1:length(pc_width)
            frac=pc_loc(k)-floor(pc_width(k)/2):floor(pc_width(k)/2)+pc_loc(k);
            frac=mod(frac-1,bins)+1; %circular wrap around indexing thinggy
            out_field=[out_field frac];
        end
        out_field=setxor(1:bins, out_field);
        out_field(isnan(out_field))=[];
    end

idx=isnan(pc_loc);
pc_loc(idx)=[];
pc_width(idx)=[];

if(~isempty(pc_width))
    reject='';
end

end