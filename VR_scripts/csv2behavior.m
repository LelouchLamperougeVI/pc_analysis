function behavior = csv2behavior(fn)
% Convert csv obtained from HaoRan's modified tabletVR to HaoRan's behavior
% format
% TODO: Should make it compatible with Dun's behavior format, which would
% require synchronization with the abf file

if nargin < 1
    C = parseCSV;
else
    C = parseCSV(fn);
end

session_cmp = @(ss, sf) strcmp(ss{2},sf{2}) & strcmp(ss{3},sf{3}) & ( str2double(ss{1}) < str2double(sf{1}) );

idx = diff(C.virtualposition(:,2));
idx = idx < -max(abs(idx))/2;
idx = find(idx);
pos_cum = C.virtualposition(:,2);
for i = 1:length(idx)
    pos_cum(idx(i)+1:end) = pos_cum(idx(i)+1:end) - pos_cum(idx(i)+1) + pos_cum(idx(i));
end
vel = diff(pos_cum)./diff(C.virtualposition(:,1));
vel = [vel(1); vel];

for i = 1:size(C.sessionStart,1)
    disp('Currently extracting:');
    disp([num2str(i) ': Animal: ' C.sessionStart{i,2} ' Session: ' C.sessionStart{i,3}]);
    
    count = 1;
    while count <= size(C.sessionFinish,1) && ~session_cmp(C.sessionStart(i,:), C.sessionFinish(count,:))
        count = count+1;
    end
    if count < size(C.sessionFinish,1)
        warning('Did not find a matching ''sessionFinish'' flag for the current session. This session will be skipped.');
        continue
    end
    
    behavior{i}.frame_ts = C.frame(C.frame(:,1) >= str2double(C.sessionStart{i,1}) & C.frame(:,1) <= str2double(C.sessionFinish{count,1}),1)';
    behavior{i}.trials_ts = knnsearch(behavior{i}.frame_ts', C.reward);
    behavior{i}.trials = behavior{i}.frame_ts(behavior{i}.trials_ts);
    idx = knnsearch(C.virtualposition(:,1), behavior{i}.frame_ts');
    behavior{i}.unit_pos = C.virtualposition(idx,2)';
    behavior{i}.unit_vel = vel(idx)';
end

