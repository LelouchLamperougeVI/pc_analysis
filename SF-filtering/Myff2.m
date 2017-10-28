function [Y fx fy mfx mfy] = Myff2(X, m, n)
% perform fft2 & fftshift, and calculate frequency vectors and matrices
% X : data to perform fft2
% m : number of rows to perform fft2
% n : number of columns to perform fft2
% Y : result of fft2 and fftshift operation for X
% fx : frequency vector
% fy : frequency vector
% mfx : frequency matrix
% mfy : frequency matrix

if nargin == 1
    [m n] = size(X);
end

% 2d fft
Y = fft2(X, m, n);
Y = fftshift(Y);

% obtain frequency (cycles/pixel)
f0 = floor([m n] / 2) + 1;

% notes : imagesc (matlab function to display image) with axis xy
%   horizontal (x) axis (low -> high) : matrix column (left -> right)
%   vertical (y) axis (low -> high) : matrix row (bottom -> top)
fy = ((m: -1: 1) - f0(1) + 1) / m;
fx = ((1: n) - f0(2)) / n;

[mfx mfy] = meshgrid(fx, fy);
