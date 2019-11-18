function deconv=burst_shuffler(deconv,fs,thres)
% The rational behind the standard circular shuffling procedure is to
% disrupt the relationship between spikes and animal behavior, while
% preserving the temporal characteristics within the spike trains.
%
% Problem arises when the animal behavior is static periodic.
%
% The proposed method preserves ISI for burst events, but shuffle
% everything else randomly.
%
% Inputs: 
%   deconv (assume filtered)
%   fs: sampling frequency
%   thres: single/burst threshold (default 1s ~fall time of GCaMP6s)
%
% Update (2019-11-14): this shit needed some insane optimization; the code
% is really convoluted now... to the future me that identified a bug: I'm
% sorry. Check version history to gain an understanding of the original
% concept.

if nargin<2
    fs=19;
end
if nargin<3
    thres=1;
end
thres=fs*thres;

for i=1:size(deconv,2)
    spk=find(deconv(:,i)>0);
    if isempty(spk)
        continue;
    end
    ISI=diff(spk);
    single=ISI>thres;
    single=find([true;single]);
    
    burst_ind = [ spk(single), spk([single(2:end) - 1; end]) ];
    burst_frames = size(deconv,1) - sum( diff(burst_ind, 1, 2) + 1 ); % number of frames to pad
    pad = rand(length(single) + 1,1);
    pad = pad ./ sum(pad) .* (burst_frames - (length(single) + 1) * thres) + thres;
    pad = round(pad);
    pad(end-1) = pad(end-1) + burst_frames - sum(pad);
    pad(end) = [];
    
    rorder = randperm(length(single));
    single = [single; length(spk)+1];
    
    rspk = zeros(length(spk), 1);
    rspk(single(rorder(1)) : single(rorder(1) + 1) - 1) = spk(single(rorder(1)) : single(rorder(1) + 1) - 1) - spk(single(rorder(1))) + pad(1) + 1;
    for ii = 2:length(single) - 1
        rspk(single(rorder(ii)) : single(rorder(ii) + 1) - 1) = ...
            spk(single(rorder(ii)) : single(rorder(ii) + 1) - 1) - spk(single(rorder(ii))) + pad(ii) + rspk(single(rorder(ii-1) + 1) - 1) + 1;
    end
    
    temp = zeros(size(deconv,1),1);
    temp(rspk) = deconv(spk,i);
    deconv(:,i) = temp;
end