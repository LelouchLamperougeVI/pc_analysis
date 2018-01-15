function lopes_pca(deconv,smooth)
% PCA cell-assemblies detection methods adapted from Lopes-dos-Santos et
% al. 2011
deconv=fast_smooth(deconv,smooth);
deconv=(deconv-mean(deconv))./std(deconv); % zscore each neuron's spike activity
corrM=corr(deconv); % compute autocorrelation matrix
[V,e]=eig(corrM); % get eigenvalues and vectors/PCs

thres=MPdist(size(deconv,1),size(deconv,2)); % get upper/lower thresholds for eigenvalue significance
e=sum(e); % vectorize eigenvalues matrix
assemblies=e>thres(2); % identify PC's associated with assemblies
numNeurons=sum(e<thres(1) | e>thres(2)); % determine number of neurons participating in assemblies
assembly_space=V(:,assemblies); % PCs with significant eigenvalues
interaction=assembly_space*assembly_space'; % interaction matrix
vect_lengths=arrayfun(@(x) norm(assembly_space(x,:)),1:size(assembly_space,1));
[~,idx]=sort(vect_lengths);
assembly_neurons=idx(end-numNeurons+1:end); % neurons with biggest norms in assembly space
interaction=interaction(assembly_neurons,:); % remove non-assembly neurons from interaction matrix
interaction=interaction(:,assembly_neurons);
diagonal=repmat(vect_lengths(assembly_neurons)',1,size(interaction,1));
interaction=interaction./diagonal'; % final interaction matrix is the scalar projection of neurons onto each other in assembly space