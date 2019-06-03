function detrend(obj,wdw)
% detrend signal using window of length wdw at half steps (default 1 sec)
if nargin < 2
    wdw=1;
end
obj.lfp = locdetrend(obj.lfp.lfp,obj.lfp.fs,[wdw wdw/2]);