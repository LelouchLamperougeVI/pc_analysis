% Load CSV file
behavior=categories_import;

% Manually load deconv

% Visualize CSV
plot_categories(behavior,deconv);

% Plot for individual neurons
analysis=cat_analysis(behavior,deconv);
plot_cat_analysis(analysis);

% Filter out left firing rate component
signal=ca_filt(deconv);

% Binarize events
signal=double(signal>0);

% Construct q-matrix
bins=10;

idx=1:bins:size(signal,1);
qMatrix=arrayfun(@(x) sum(signal(idx(x):idx(x+1),:))', 1:length(idx)-1,'uniformoutput',false);
qMatrix=cell2mat(qMatrix)';

% Run Lopes-dos-Santos PCA
[assemblies,R]=lopes_pca(qMatrix,0,1);
