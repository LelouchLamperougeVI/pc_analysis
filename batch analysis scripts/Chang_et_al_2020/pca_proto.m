clear all

% test = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');
test = ensemble('/mnt/storage/rrr_magnum/M2/EE006/2019_05_30/2019_05_30_3.abf');
test.set_ops('e_size',5);
test.set_ops('clust_method','thres');
test.set_ops('sig', .2);
test.remove_mvt;
test.cluster;

test.hPICA;
test.plot('hPICA');

