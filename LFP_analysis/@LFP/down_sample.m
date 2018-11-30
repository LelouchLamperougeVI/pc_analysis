function down_sample(obj,Fs)
% Downsample

factor = floor(obj.fs/Fs);
obj.fs = obj.fs * ceil(length(obj.lfp)/factor)/length(obj.lfp);
obj.lfp = decimate(obj.lfp,factor,'fir');

end