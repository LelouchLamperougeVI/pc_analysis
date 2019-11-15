function down_sample(obj,Fs)
% Downsample

if isempty(obj.lfp.lfp)
    return
end

factor = floor(obj.lfp.fs/Fs);
if ~factor
    error('requested sampling rate is too fast');
end
obj.lfp.lfp = decimate(obj.lfp.lfp,factor,'fir');
obj.lfp.fs = obj.lfp.fs * ceil(length(obj.lfp.lfp)/factor)/length(obj.lfp.lfp);
obj.lfp.ts = linspace(0,(length(obj.lfp.lfp) - 1)/obj.lfp.fs,length(obj.lfp.lfp));

end