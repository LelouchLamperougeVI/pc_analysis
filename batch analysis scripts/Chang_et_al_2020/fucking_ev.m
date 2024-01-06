clear all
load('/mnt/storage/rrr_magnum/M2/cross_days.mat', 'EV')
ev1 = EV;
load('/mnt/storage/HaoRan/RRR_motor/M2/cross_days.mat', 'EV')
ev2 = EV;
EV = [ev2; ev1];
signrank(diff(EV, 1, 2), 0, 'tail', 'left')

figure
boxplot(EV)