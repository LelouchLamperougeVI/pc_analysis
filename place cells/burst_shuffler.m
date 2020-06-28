function [ret_deconv, spk] = burst_shuffler(deconv,fs,thres, spk)
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
%   spk: reuse spk in shuffle test for improved performance
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

ret_deconv = zeros(size(deconv));

if nargin < 4
    spk = cell(size(deconv,2), 1);
end

for i=1:size(deconv,2)
    if nargin < 4
        spk{i} = find(deconv(:,i)>0);
    end
    if isempty(spk{i})
        continue;
    end
    ISI=diff(spk{i});
    single=ISI>thres;
    single=find([true;single]);
    
    burst_ind = [ spk{i}(single), spk{i}([single(2:end) - 1; end]) ];
    burst_frames = size(deconv,1) - sum( diff(burst_ind, 1, 2) + 1 ); % number of frames to pad
    burst_spks = [diff(single); length(spk{i}) - single(end) + 1]; % number of spks per burst
    pad = rand(length(single) + 1,1);
    pad = pad ./ sum(pad) .* (burst_frames - (length(single) + 1) * thres) + thres;
    pad = round(pad);
    pad(end-1) = pad(end-1) + burst_frames - sum(pad);
    pad(end) = [];
    
    rorder = randperm(length(single));
    [~, reverse] = sort(rorder);

    rburst_ind = burst_ind(rorder, :);
    rburst_len = diff(rburst_ind, 1, 2);
    ind = cumsum( [-1; rburst_len(1:end-1)] + 1 ) + cumsum(pad);

    shift = repelem(ind(reverse) - burst_ind(:, 1), burst_spks);
    
    ret_deconv(spk{i} + shift, i) = deconv(spk{i}, i);
end