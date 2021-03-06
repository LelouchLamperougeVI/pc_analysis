function [registered, stats] = regmasks(fixed, moving, varargin)
% Custom rigid registration algorithm for ROI masks.
% Based on spatial cross-correlation with mini-batch gradient descent
% optimizer. Mini-batch was mainly used to prevent converging onto local
% maxima...
%
% Usage: registered = regmasks(fixed, moving, varargin)
%
% positional arguments:
%   fixed               fixed reference mask
%   moving              moving mask
%
% optional arguments:
%   lag, 10             maximum lag for spatial cross-correlation
%   theta, 1            initial guesstimate of theta angle (in degrees)
%   gamma, 1            learning rate
%   max_it, 1e3         maximum number of iterations
%   precision, 1e-9     fault tolerance
%   percentile, 15      percentage of moving ROIs to use
%   all (false)         use all ROIs
%   verbose (false)     for debugging and performance optimization
%
% output:
%   registered          registered moving mask

ops = parse_args(varargin);

amat = @(theta) [ cosd(theta)   sind(theta)     0; % obtain affine transformation matrix given angle (degrees)
                 -sind(theta)   cosd(theta)     0
                            0             0     1 ];
%% step 1: find ROI centroids %%
f_centroids = zeros(max(fixed(:)), 2); % centroids of ROIs in fixed plane
for i = 1 : max(fixed(:))
    props = regionprops(fixed==i);
    f_centroids(i,:) = props.Centroid;
end

m_centroids = zeros(max(moving(:)), 2); % centroids of ROIs in moving plane
for i = 1 : max(moving(:))
    props = regionprops(moving==i);
    m_centroids(i,:) = props.Centroid;
end

d = m_centroids - permute(f_centroids, [3 2 1]);
d = sqrt(sum(d.^2, 2));
d = min(d, [], 3);

if ops.all; ops.percentile = 100; end
neurIDs = find( d < prctile(d, ops.percentile) )';

%% step 2: mini-batch gradient descent %%
% turns out mini-batch was too unstable...
mask = double(ismember(moving, neurIDs));
fixed = double(~~fixed); % convert template to binary image

temp = bxcorr2(fixed, mask, ops.lag, 'unbiased');
r_max = max(temp(:)) / sum(mask(:)) * numel(mask);

if ops.verbose; figure; hold on; end

step = ops.theta;
stats.theta = ops.theta;
stats.it = 0;
while stats.it < ops.max_it
    tform = affine2d( amat(stats.theta) );
    transformed = imwarp(mask, tform, 'nearest');
    coor = ( size(transformed) - size(fixed) ) ./ 2;
%         transformed = transformed(floor(coor(1)) : end - ceil(coor(1)) - 1, floor(coor(1)) : end - ceil(coor(1)) - 1);
    transformed = transformed( floor(coor(1)) + 1 : end - ceil(coor(1)), floor(coor(2)) + 1 : end - ceil(coor(2)));
    
    temp = bxcorr2(fixed, transformed, ops.lag, 'unbiased');
    r = max(temp(:)) / sum(transformed(:)) * numel(transformed);
    step = ops.gamma * ( r - r_max ) / step;
    stats.theta = stats.theta + step;
    r_max = r;
    
    stats.it = stats.it + 1;
    if abs(step) < ops.precision
        break;
    end
    
    if ops.verbose; plot3(stats.it, r, stats.theta, 'ko'); disp([ 'next theta: ' num2str(stats.theta) ' r_max:' num2str(r_max) ' step:' num2str(step) ]); end
end

if ops.verbose; xlabel('iteration'); ylabel('r_{max} (mask \mu overlap)'); zlabel('\theta'); title(['learning rate \gamma = ' num2str(ops.gamma)]); view(3); end

%% step 3: final translation transform %%
tform = affine2d( amat(stats.theta) );
registered = imwarp(moving, tform, 'nearest');
coor = ( size(registered) - size(fixed) ) ./ 2;
registered = registered( floor(coor(1)) + 1 : end - ceil(coor(1)), floor(coor(2)) + 1 : end - ceil(coor(2)));
[r, lags] = bxcorr2(fixed, registered, 30, 'unbiased');
stats.r = max(r(:));
[x, y] = find(r == stats.r);

if lags.x(x) > 0
    registered = padarray(registered, [lags.x(x) 0], 'pre');
    registered(end-lags.x(x)+1 : end, :) = [];
elseif lags.x(x) < 0
    registered = padarray(registered, [-lags.x(x) 0], 'post');
    registered(1 : -lags.x(x), :) = [];
end

if lags.y(y) > 0
    registered = padarray(registered, [0 lags.y(y)], 'pre');
    registered(:, end-lags.y(y)+1 : end) = [];
elseif lags.y(y) < 0
    registered = padarray(registered, [0 -lags.y(y)], 'post');
    registered(:, 1 : -lags.y(y)) = [];
end

if ops.verbose
    figure;
    subplot(1,3,1);
    imshowpair(~~fixed,~~moving)
    title('before')
    subplot(1,3,2);
    imshowpair(~~fixed,~~registered)
    title('after')
    subplot(1,3,3);
    imshow(mask)
    title('chosen')
end

function ops = parse_args(inputs)

ops.lag = 10;
ops.gamma = 1;
ops.theta = 1;
ops.max_it = 1e3;
ops.precision = 1e-9;
ops.percentile = 15;
ops.all = false;
ops.verbose = false;

count = 1;
while count <= length(inputs)
    switch lower(inputs{count})
        case {'lag', 'range', 'max_lag'}
            ops.lag = inputs{count + 1};
        case {'gamma', 'rate', 'learningrate'}
            ops.gamma = inputs{count + 1};
        case {'init', 'initial', 'theta'}
            ops.theta = inputs{count + 1};
        case {'iterations', 'max_it'}
            ops.max_it = inputs{count + 1};
        case {'precision', 'tolerance'}
            ops.precision = inputs{count + 1};
        case {'percentile', 'prctile'}
            ops.percentile = inputs{count + 1};
        case {'all'}
            ops.all = true;
            count = count - 1;
        case {'verbose', 'plot'}
            ops.verbose = true;
            count = count - 1;
        otherwise
            error(['The input ''' inputs{count} ''' is not a valid argument.']);
    end
    count = count + 2;
end