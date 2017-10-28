% spatial frequency (SF) filtering by SF-bandpass filter

% the SF filter is unselective to orientation (doughnut-shaped in the SF
% domain).

filename = 'lol.tif';

% parameters for SF filtering
ctrSF = 0.01; 0.06; % center SF
s = 0.04; % sigma of Gaussian function

show_log_amp = true; % flag to show log SF spectrum instead of SF spectrum
min_amp_to_show = 10 ^ -10; % small positive value to replace 0 for log SF spectrum

L = GetLuminanceImage(filename);

% calculate the number of points for FFT (power of 2)
FFT_pts = 2 .^ ceil(log2(size(L)));

[A fx fy mfx mfy] = Myff2(L, FFT_pts(1), FFT_pts(2));

SF = sqrt(mfx .^ 2 + mfy .^ 2);

% SF-bandpass and orientation-unselective filter
filt = exp(-(SF - ctrSF) .^ 2 / (2 * s ^ 2));

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
imagesc(fx, fy, A);
axis xy;
axis square;
set(gca, 'TickDir', 'out');
title('amplitude spectrum');
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
