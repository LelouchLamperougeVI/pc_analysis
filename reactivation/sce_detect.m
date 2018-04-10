function sce=sce_detect(deconv, sig)
% Detect the locations of SCEs using a PCA method
% Description: Get significant eigenvectors using MP distribution and
% compute the average vector normalized by eigenvalues
% Inputs:
%   deconv: deconvolved signal
%   sig: standard deviation for gausian smoothing
% Output:
%   sce: averaged eigenvector

signal=fast_smooth(deconv,sig);
signal=zscore(signal);
M=cov(signal);
[e,v]=eig(M);