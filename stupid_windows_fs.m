function deconv=stupid_windows_fs(deconv)

frames_per_file=2000;

numFiles=size(deconv,1)/frames_per_file;

excess=numFiles-9;
excess=excess*frames_per_file + frames_per_file;

deconv=[deconv(1:2000,:); deconv(excess+1:end,:); deconv(2001:excess,:)];