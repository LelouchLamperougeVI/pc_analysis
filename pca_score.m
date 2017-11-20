function scores=pca_score(deconv,sig)
% get PCA scores for resting state series

d=fast_smooth(deconv,19.1*sig);
d=zscore(d);
[~,scores]=pca(d);