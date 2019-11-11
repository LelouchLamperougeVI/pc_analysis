clear all
% 
maskA = load('/mnt/storage/data/EE004/2019_05_01/1/plane1/masks_neurons.mat');
maskA = double(~~maskA.maskNeurons);
maskB = load('/mnt/storage/data/EE004/2019_05_06/1/plane1/masks_neurons.mat');
maskB = double(~~maskB.maskNeurons);
figure;
imshowpair(maskA,maskB)


%%
lag = 30;
rot_range = 5;
step = 0.1;

lol = [];

for theta = -rot_range:step:rot_range
    tform = [ cosd(theta)   sind(theta)     0;
             -sind(theta)   cosd(theta)     0
                        0             0     1 ];
    tform = affine2d(tform);
    transformed = imwarp(maskB, tform, 'nearest');
    coor = ( size(transformed) - size(maskA) ) ./ 2;
    try
        transformed = transformed(floor(coor(1)) : end - ceil(coor(1)) - 1, floor(coor(1)) : end - ceil(coor(1)) - 1);
    catch
    end
    
    r = bxcorr2(maskA, transformed, lag, 'unbiased');
    lol = [lol; theta, max(r(:))]
end

%%
theta = lol(lol(:,2)==max(lol(:,2)), 1);
tform = [ cosd(theta)   sind(theta)     0;
         -sind(theta)   cosd(theta)     0
                    0             0     1 ];
tform = affine2d(tform);
transformed = imwarp(maskB, tform, 'nearest');
coor = ( size(transformed) - size(maskA) ) ./ 2;
try
    transformed = transformed(floor(coor(1)) : end - ceil(coor(1)) - 1, floor(coor(1)) : end - ceil(coor(1)) - 1);
catch
end

[r, lags] = bxcorr2(maskA, transformed, 30, 'unbiased');
m = max(r(:));
[x, y] = find(r == m);
lags.x(x)
lags.x(y)

figure;
imshowpair(maskA,maskB)

if lags.x(x) > 0
    transformed = padarray(transformed, [lags.x(x) 0], 'pre');
    transformed(end-lags.x(x)+1 : end, :) = [];
elseif lags.x(x) < 0
    transformed = padarray(transformed, [-lags.x(x) 0], 'post');
    transformed(1 : -lags.x(x), :) = [];
end

if lags.y(y) > 0
    transformed = padarray(transformed, [0 lags.y(y)], 'pre');
    transformed(:, end-lags.y(y)+1 : end) = [];
elseif lags.y(y) < 0
    transformed = padarray(transformed, [0 -lags.y(y)], 'post');
    transformed(:, 1 : -lags.y(y)) = [];
end

figure;
imshowpair(maskA,transformed)

bxcorr2(maskA, transformed, 3, 'unbiased')

