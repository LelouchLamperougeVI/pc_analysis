function s = strrank(s)
% Sometimes, strings are unordered.
% This could be caused by the stupid filesystem
% E.g. plane1 plane10 plane11 plane2 plane3 ...
% This function takes a cell array of strings and order by the last
% numerical value

idx = regexp(s,'\d+', 'match');
idx = cellfun(@(x) str2double(x{end}), idx);
[~,idx] = sort(idx);

s = s(idx);