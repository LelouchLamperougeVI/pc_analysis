function I = fast_smooth2(A, varargin)
% 2D version of fast_smooth with more sophisticated methodology (as I'm
% getting better at programming matlab, or so I believe)
% As with before, accounts for edge underestimation and NaNs
% Inputs:
%   A: matrix to convolve
% Name, Value:
%   'method':   'gauss' (default) or 'mean'
%   'sig':      1 (default)
%   'wdw':      5 (default)

gausswin2 = @(x, y, sig) exp( - (x.^2 + y.^2) ./2 ./(sig.^2) );

ops = parse_inputs(varargin);

switch ops.method
    case 'gauss'
        kernelSize = ceil(10 * ops.sig);
        kernelSize = kernelSize - ~mod(kernelSize,2); %avoid even kernel sizes
        
        x = repmat(-floor(kernelSize/2) : floor(kernelSize/2), kernelSize, 1);
        y = x';
        
        kernel = gausswin2(x, y, ops.sig);
        kernel = kernel ./ sum(kernel(:));
    case 'mean'
        kernelSize = ops.wdw;
        kernel = ones(kernelSize) ./ ( kernelSize .^2 );
end

nan_idx = isnan(A);
A(nan_idx) = 0;
taper = padarray(~nan_idx, [floor(kernelSize/2) floor(kernelSize/2)], 0, 'both');
A = padarray(A, [floor(kernelSize/2) floor(kernelSize/2)], 0, 'both');

I = conv2(A, kernel, 'same') ./ conv2(taper, ones(kernelSize)./(kernelSize.^2), 'same'); % divide by convolution using unit kernel
I([1:floor(kernelSize/2) end-floor(kernelSize/2)+1:end],:) = [];
I(:, [1:floor(kernelSize/2) end-floor(kernelSize/2)+1:end]) = [];

I(I==0) = nan;


function ops = parse_inputs(in)
% name, value pairs parser
ops.method = 'gauss'; % choice between 'gauss' Gaussian kernel and 'mean' averaging kernel
ops.sig = 1; % sigma for gaussian kernel
ops.wdw = 5; % size of averaging kernel

idx = 1;
while idx < length(in)
    switch lower(in{idx})
        case 'method'
            ops.method = in{idx+1};
        case 'sig'
            ops.sig = in{idx+1};
        case 'wdw'
            ops.wdw = in{idx+1};
        otherwise
            warning(['''' in{idx} ''' is not a valid parameter. Ignored']);
    end
    idx = idx+2;
end