function registered = regmasks(fixed, moving)
% Custom rigid registration algorithm for ROI masks
% Based on spatial cross-correlation with gradient descent optimizer

lag = 30;
gamma = 1;
theta = 0.01;
max_it = 100;
precision = 1e-9;

amat = @(theta) [ cosd(theta)   sind(theta)     0; % obtain affine transformation matrix given angle (degrees)
                 -sind(theta)   cosd(theta)     0
                            0             0     1 ];

alpha = min([sum(fixed(:)) sum(moving(:))]) / numel(fixed);
r = bxcorr2(fixed, moving, lag, 'unbiased');
r_max = max(r(:));
it = 0;
while it < max_it
    tform = affine2d( amat(theta) );
    transformed = imwarp(fixed, tform, 'nearest');
    coor = ( size(transformed) - size(fixed) ) ./ 2;
    try
        transformed = transformed(floor(coor(1)) : end - ceil(coor(1)) - 1, floor(coor(1)) : end - ceil(coor(1)) - 1);
    catch
    end
    
    r = bxcorr2(fixed, transformed, lag, 'unbiased');
    r = max(r(:));
    step = gamma * ( r - r_max ) / alpha;
    r_max = r;
    theta = theta + step;
    
    it = it + 1;
    disp([ 'theta: ' num2str(theta) ' r_max:' num2str(r_max) ' step:' num2str(step) ]);
    if abs(step) < precision
        break;
    end
end