precision=5;

assemblies=cluster_mi(deconv,'prune',15,'plotFlag',true,'shuffle',100,'sig',5,'precision',precision);

[sce,loss]=mi_sce(assemblies,deconv,precision);

plot_assemblies(assemblies,sce,deconv);