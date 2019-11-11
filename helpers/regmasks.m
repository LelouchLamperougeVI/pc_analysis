function registered = regmasks(fixed, moving)
% Custom rigid registration algorithm for ROI masks
% Based on spatial cross-correlation with gradient descent optimizer
% Canonical GD is not suitable for correlation values (some weird-ass
% convoluted methods are used instead...)

lag = 30;
gamma = 1;
theta = 0.1;
max_it = 20;
precision = 1e-9;
min_step = 0.01;

amat = @(theta) [ cosd(theta)   sind(theta)     0; % obtain affine transformation matrix given angle (degrees)
                 -sind(theta)   cosd(theta)     0
                            0             0     1 ];

alpha = min([sum(fixed(:)) sum(moving(:))]) / numel(fixed);
r = bxcorr2(fixed, moving, lag, 'unbiased');
r_max = max(r(:));
it = 0;
while it < max_it
    tform = affine2d( amat(theta) );
    transformed = imwarp(moving, tform, 'nearest');
    coor = ( size(transformed) - size(fixed) ) ./ 2;
    try
        transformed = transformed(floor(coor(1)) : end - ceil(coor(1)) - 1, floor(coor(1)) : end - ceil(coor(1)) - 1);
    catch
    end
    
    r = bxcorr2(fixed, transformed, lag, 'unbiased');
    r = max(r(:));
    ( r - r_max ) / alpha
    step = gamma * ( r - r_max ) / alpha;
    r_max = r;
    theta = theta + max([abs(step) min_step]) * sign(step);
    
    it = it + 1;
    disp([ 'next theta: ' num2str(theta) ' r_max:' num2str(r_max) ' step:' num2str(step) ]);
    if abs(step) < precision
        break;
    end
end

[r, lags] = bxcorr2(fixed, transformed, 30, 'unbiased');
m = max(r(:));
[x, y] = find(r == m);

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
subplot(1,2,1);
title('before')
imshowpair(fixed,moving)
subplot(1,2,2);
title('after')
imshowpair(fixed,transformed)
