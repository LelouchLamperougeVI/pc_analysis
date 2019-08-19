% Compile all source codes on Linux system

mex -g -O mat_conv.c convolve.c
mex -g -O faster_smooth.c convolve.c
mex -g -O mat_circshift.c
