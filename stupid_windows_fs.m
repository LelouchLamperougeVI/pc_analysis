function deconv=stupid_windows_fs(deconv)
% somewhere along the preprocessing, the stupid file system decided to
% order file numbers in the following way: 1, 10, 11, 12, ... 2, 3, 4, ...
% as a result, the deconv is all messed up
% pass deconv through this function to fix everything :)

frames_per_file=2000;

numFiles=size(deconv,1)/frames_per_file;

excess=numFiles-9;
excess=excess*frames_per_file + frames_per_file;

deconv=[deconv(1:2000,:); deconv(excess+1:end,:); deconv(2001:excess,:)];