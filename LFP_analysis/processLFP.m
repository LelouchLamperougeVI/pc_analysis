function lfp=processLFP(signal,fs)
% Detrend, denoise and filter lfp signal
%
% HaoRan Chang


filt60 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);

filtSWR=designfilt('bandpassfir','stopbandfrequency1',149,'passbandfrequency1',150,'passbandfrequency2',250,'stopbandfrequency2',251,...
                'stopbandattenuation1',40,'passbandripple',1,'stopbandattenuation2',30,...
                'designmethod','kaiserwin','samplerate',fs);

signal=detrend(signal,'linear',0.05*fs);
signal=filtfilt(filt60,signal);
signal=filtfilt(filtSWR,signal);