function assemblies=bin_pca(q)
  % PCA cell-assembly detection method for Q-matrix

  q=zscore(q); %normalize trains

  corrM=corr(q); % compute autocorrelation matrix
  [V,e]=eig(corrM); % get eigenvalues and vectors/PCs
  e=diag(e); % vectorize eigenvalues matrix

  [thres,epdf]=MPdist(size(q,1),size(q,2),min(e):0.01:max(e)); % get upper/lower thresholds for eigenvalue significance

  assemblies.accuracy=sum(e<thres(2) & e>thres(1))/length(e);
