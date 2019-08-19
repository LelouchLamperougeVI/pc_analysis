function deconv=burst_shuffler(deconv,fs,thres)
% The rational behind the standard circular shuffling procedure is to
% disrupt the relationship between spikes and animal behavior, while
% preserving the temporal characteristics within the spike trains.
%
% Problem arises when the animal behavior is static periodic.
%
% The proposed method preserves ISI for burst events, but shuffle
% everything else randomly.
% Pseudo-code:
%   get histogram of ISI between
%
% Inputs: 
%   deconv (assume filtered)
%   fs: sampling frequency
%   thres: single/burst threshold (default 1s ~fall time of GCaMP6s)
if nargin<2
    fs=19;
end
if nargin<3
    thres=1;
end
thres=fs*thres;

for i=1:size(deconv,2)
    idx=find(deconv(:,i)>0);
    if isempty(idx)
        continue;
    end
    ISI=diff(idx);
    single=ISI>thres;
    single=find([true;single]);
    
    dt=max([idx(single(2:end)-1)-idx(single(1:end-1)); idx(end)-idx(single(end))])+1;
    if floor(size(deconv,1)/dt) > length(single)
        samples=randsample(floor(size(deconv,1)/dt),length(single));
    else
        dt=1;
        single=1:length(idx);
        samples=randsample(size(deconv,1),length(single));
    end
        
    temp=zeros(size(deconv,1),1);
    for j=1:length(samples)-1
        ind=idx(single(j)):idx(single(j+1)-1);
        temp(((samples(j)-1)*dt+1):((samples(j)-1)*dt+1)+length(ind)-1)=deconv(ind,i);
    end
    ind=idx(single(end)):idx(end);
    temp(((samples(end)-1)*dt+1):((samples(end)-1)*dt+1)+length(ind)-1)=deconv(ind,i);
    deconv(:,i)=temp;
end