function down_sample(obj,Fs)
% Downsample

factor = floor(obj.fs/Fs);
if ~factor
    error('requested sampling rate is too fast');
end
obj.lfp = decimate(obj.lfp,factor,'fir');
obj.fs = obj.fs * ceil(length(obj.lfp)/factor)/length(obj.lfp);
obj.t = linspace(0,(length(obj.lfp) - 1)/obj.fs,length(obj.lfp));

end