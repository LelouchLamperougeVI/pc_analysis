function lfp=processLFP(signal,fs,d_fs)
% Downsample data and lowpass filter @ <300 Hz
%
if nargin<3
    d_fs=1250;
end

filt60 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);

filtSWR=designfilt('bandpassfir','stopbandfrequency1',149,'passbandfrequency1',150,'passbandfrequency2',250,'stopbandfrequency2',251,...
                'stopbandattenuation1',40,'passbandripple',1,'stopbandattenuation2',30,...
                'designmethod','kaiserwin','samplerate',fs);

signal=detrend(signal,'linear',0.05*fs);
signal=filtfilt(filt60,signal);
signal=filtfilt(filtSWR,signal);

