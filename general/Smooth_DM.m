function smoothed = Smooth(data,smooth)

%Smooth - Smooth using a Gaussian kernel.
%
%  USAGE
%
%    smoothed = Smooth(data,smooth)
%
%    data           data to smooth
%    smooth         vertical and horizontal standard deviations [Sv Sh]
%                   for Gaussian kernel, measured in number of samples
%                   (0 = no smoothing)

% Copyright (C) 2004-2011 by Michael Zugaro
% Modified by Dun Mao
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Vectors must be 'vertical'
if size(data,1) == 1,
	data = data';
end

if nargin < 2,
	error('Incorrect number of parameters');
end

% If Sh = Sv = 0, no smoothing required
if all(smooth==0),
	smoothed = data;
	return
end

if length(smooth) == 1,
	% For 2D data, providing only one value S for the std is interpreted as Sh = Sv = S
	smooth = [smooth smooth];
end

% Build Gaussian kernels
[vSize,~] = size(data);
% Vertical kernel
vKernelSize = min([vSize 1001]);
r = (-vKernelSize:vKernelSize)'/vKernelSize;
vKernelStdev = smooth(1)/vKernelSize;
vKernel = exp(-r.^2/(vKernelStdev+eps)^2/2);
vKernel = vKernel/sum(vKernel);

% Vector smoothing
% Prepend/append data to limit edge effects
top = flipud(data(1:vKernelSize));
bottom = flipud(data(end-vKernelSize+1:end));
data = [top;data;bottom];
% Convolve (and return central part)
tmp = conv(vKernel,data);
n = size(tmp,1);
d = n - vSize;
start = d/2+1;
stop = start + vSize - 1;
smoothed = tmp(start:stop,:);
