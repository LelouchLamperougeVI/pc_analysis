function convolved=fast_conv(x,y)
% faster convolution based on convolution theorem (only works with large
% data)

n = length(x) + length(y) - 1;
convolved = ifft(fft(x,n) .* fft(y,n));