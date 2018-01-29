function convolved=fast_conv(x,y)
% faster convolution based on convolution theorem (only works with large
% data)

n=length(x)+length(y)-1;
a=fft(x,n);
b=fft(y,n);
convolved=a.*b;
convolved=ifft(convolved);