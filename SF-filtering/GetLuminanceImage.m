function L = GetLuminanceImage(filename)
% obtain gray scale image in file
% filename : filename of image
% L : gray scale image

I = find(filename == '.', 1, 'last');
fmt = filename(I + 1: end);

[X map] = imread(filename, fmt);

X = double(X);

if ndims(X) == 3 % color image
    % obtain luminance from RGB channel values
    L = 0.30 * X(:, :, 1) + 0.59 * X(:, :, 2) + 0.11 * X(:, :, 3);
else % monochrome image
    L = X;
end

L = L - 127.5;
