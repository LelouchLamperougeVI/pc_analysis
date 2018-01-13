function lopes_pca(deconv)
% PCA cell-assemblies detection methods adapted from Lopes-dos-Santos et
% al. 2011

deconv=(deconv-mean(deconv))./std(deconv); % zscore each neuron's spike activity
corrM=corr(deconv); % compute autocorrelation matrix
[V,e]=eig(corrM); % get eigenvalues and vectors/PCs