% spatial frequency (SF) filtering by SF- and orientation-selective filter

% this filter is composed of two Gaussian functions that are symmetric
% about the origin in the SF domain. this corresponds to Fourier transform 
% of a 2D Gabor function in the space domain. therefore, filtering operation 
% by this filter is equivalent to computation performed in simple cells in 
% the early visual cortex.

filename = 'lol.tif';

% parameters for SF filtering
ctrSF = 0.06 * [cos(pi / 6) sin(pi / 6)]; % center SF
s = 0.02; % sigma of Gaussian function

show_log_amp = true; % flag to show log SF spectrum instead of SF spectrum
min_amp_to_show = 10 ^ -10; % small positive value to replace 0 for log SF spectrum 

L = GetLuminanceImage(filename);

% calculate the number of points for FFT (power of 2)
FFT_pts = 2 .^ ceil(log2(size(L)));

[A fx fy mfx mfy] = Myff2(L, FFT_pts(1), FFT_pts(2));

% SF filter 1 : center SF (ctrSF(1), ctrSF(2))
filt1 = exp(-((mfx - ctrSF(1)) .^ 2 + (mfy - ctrSF(2)) .^ 2) / (2 * s ^ 2));

% SF filter 2 : center SF (ctrSF(1), -SFpeak(2))
filt2 = exp(-((mfx + ctrSF(1)) .^ 2 + (mfy + ctrSF(2)) .^ 2) / (2 * s ^ 2));

% SF- and orientation-selective filter
% add the above two filters to make a filter symmetric about the origin
filt = filt1 + filt2;

A_filtered = filt .* A; % SF filtering
L_filtered = real(ifft2(ifftshift(A_filtered))); % IFFT
L_filtered = L_filtered(1: size(L, 1), 1: size(L, 2));


figure(1);
clf reset;

colormap gray;

% luminance image
subplot(2, 2, 1);
imagesc(L);
colorbar;
axis square;
set(gca, 'TickDir', 'out');
title('original image');
xlabel('x');
ylabel('y');

A = abs(A);
if show_log_amp
    A(find(A < min_amp_to_show)) = min_amp_to_show; % avoid taking log 0
    A = log10(A);
end

% spectral amplitude
subplot(2, 2, 2);
imagesc(fx, fy, abs(A));
axis xy;
axis square;
title('amplitude spectrum');
set(gca, 'TickDir', 'out');
xlabel('fx (cyc/pix)');
ylabel('fy (cyc/pix)');

% filter in the SF domain
subplot(2, 2, 3);
imagesc(fx, fy, filt);
axis xy;
axis square;
set(gca, 'TickDir', 'out');
title('filter in the SF domain');
xlabel('fx (cyc/pix)');
ylabel('fy (cyc/pix)');

% filtered image
subplot(2, 2, 4);
imagesc(L_filtered);
colorbar;
axis square;
set(gca, 'TickDir', 'out');
title('filtered image');
xlabel('x');
ylabel('y');
