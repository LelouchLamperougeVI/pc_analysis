function SCE = basic_ensembles(varargin)
% A very basic analysis for detecting synchronous calcium events (SCE) by
% performing a sum of Q over zscored deconv
% This should serve as a temporary solution for the Philosophical
% Transactions paper, while we work on finalizing the MI method
% Inputs:
%   'deconv'
%   'ts':           timestamps for 2p frames in seconds
%
%   'sig':          gaussian smoothing kernel SD (default - 50 ms)
%   'thres':        threshold for SCE expressed in SDs (default - 3 SD)
%   'off_thres':    threshold for SCE onset/offset (default - 1 SD)
%   'gaps':         distance in seconds for 2 SCEs to be considered as a singular event (default - 250 ms)
%
% Output: structure 'SCE.
%   dur':       duration of a SCE event
%   on':        onset time of SCEs
%   peaks':     time of SCE response peaks

deconv=varargin{1};
ts=varargin{2};
[sig,thres,off_thres,gaps]=parse_input(varargin);

fs=1/median(diff(ts));
sig=fs*sig;
gaps=round(fs*gaps);


deconv=zscore(deconv);
deconv=fast_smooth(deconv,sig);

mua=sum(deconv,2);
thres=mean(mua)+thres*std(mua);
sce=mua>thres;
sce=fill_gaps(sce,gaps);

off_thres=mean(mua)+off_thres*std(mua);
off=mua<off_thres;
off=off.*~sce;

heads=get_head(sce);
tails=get_head(sce(end:-1:1));
tails=tails(end:-1:1);
if heads(1)==1
    heads(1)=[];
    tails(1)=[];
end
if tails(end)==length(sce)
    heads(end)=[];
    tails(end)=[];
end
heads=find(heads);tails=find(tails);

for i=1:length(heads)
    heads(i)=knnsearch();
end


function [sig,thres,off_thres,gaps]=parse_input(inputs)
sig=.05;
thres=3;
off_thres=1;
gaps=.25;

idx=3;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 'thres'
            idx=idx+1;
            thres=inputs{idx};
        case 'gaps'
            idx=idx+1;
            gaps=inputs{idx};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx=idx+1;
end