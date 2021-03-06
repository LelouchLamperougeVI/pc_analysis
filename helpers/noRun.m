function thres=noRun(unit_vel)
% find velocity threshold for non-running epochs

[c,centres]=hist((abs(unit_vel)),length(unit_vel)/2);
[~,thres]=findpeaks(-c);
try
    thres=(centres(thres(1)));
catch
    thres=inf;
end