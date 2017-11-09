function thres=noRun(unit_vel)
% find velocity threshold for non-running epochs

[c,centres]=hist((abs(unit_vel)),length(unit_vel)/2);
thres=findpeaks(-c);
thres=(centres(thres.loc(1)));