function down_sample(obj,Fs)
% Downsample

factor = floor(obj.fs/Fs);
if ~factor
    error('requested sampling rate is too fast');
end
obj.lfp = decimate(obj.lfp,factor,'fir');
obj.fs = obj.fs * ceil(length(obj.lfp)/factor)/length(obj.lfp);

end