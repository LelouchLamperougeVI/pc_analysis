clear all

root = '/mnt/storage/rrr_magnum/M2';
animal = 'EE001';
date = '2018_12_28';
session = '2';

% chan = [1 2 3 5 NaN 6 NaN];
chan = [1 2 3 6 NaN 5 NaN];

ass = ensemble(fullfile(root, animal, date, session, [date '_' session '.abf']), 'channels', chan);