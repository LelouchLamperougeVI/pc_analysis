%% Import stim and behavioral data for HD cells
clear frame_ts

% 
importFlag={'"rsc32 t1"'};
fn='VR20170722140623.csv';

importFlag={'"rsc32 t2"'};
fn='VR20170722144334.csv';
% 
% importFlag={'"rsc32"','"deep"'};
% fn='VR20170721135307.csv';
% 
% importFlag={'"rsc34"'};
% fn='VR20170710161414.csv';

% importFlag={'"rsc032 randT"'};
% fn='VR20170906171140.csv';

% importFlag={'"rsc032 contV dark"'};
% fn='VR20170906174039.csv';
% 
importFlag={'"rsc34 v360"','"dark"'};
importFlag={'"dark"'};
fn='VR20170726143105.csv';
% 
% importFlag={'"rsc32"','"deep"'};
% fn='VR20170720150520.csv';



fid=fopen(fn);
C=textscan(fid,'%s','delimiter',',');
fclose(fid);
C=C{1,1};

recStart=cellfun(@(x) strcmp(x,importFlag{1}), C);
recStart=find(recStart,1);

if length(importFlag)==2
    recFinish=cellfun(@(x) strcmp(x,importFlag{2}), C);
    recFinish=find(recFinish,1);
end

idx=cellfun(@(x) strcmp(x, 'rotation'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
for i=1:length(idx)
    count=C{idx(i)+3};
    rotation(i)=str2double(count);
    count=C{idx(i)-1};
    ts(i)=str2double(count);
end

% get frame timestamps; sync this shit with ts then interpolate
idx=cellfun(@(x) strcmp(x, 'count-21'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
for i=1:length(idx)
    count=C{idx(i)-1};
    frame_ts(i)=str2double(count);
end
frame_ts=frame_ts(1:2:end);

% get trial timestamps
idx=cellfun(@(x) strcmp(x, 'trigger'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
for i=1:length(idx)
    count=C{idx(i)-1};
    trials(i)=str2double(count);
end
trials(trials>frame_ts(end))=[];

trials_ts=arrayfun(@(x) find(frame_ts>=x,1),trials);

% export to unit_pos
unit_rot=arrayfun(@(x) find(ts>=x,1), frame_ts);
unit_rot=rotation(unit_rot);

[~,trials_ts]=findpeaks(unit_rot);

%% pin 21 correction
idx=1344;
unit_rot(1:idx)=[];
frame_ts(1:idx)=[];
idx=find(trials>frame_ts(1),1);
trials(1:idx-1)=[];

% unit_rot(end-9667:end)=[];
% frame_ts(end-9667:end)=[];
% idx=find(trials>frame_ts(end),1);
% trials(idx:end)=[];


%% find trial timestamps (independent of reward)
% run this after correcting for fake registerations on pin 21
[~,trials_ts]=findpeaks(unit_rot);
trials_ts(unit_rot(trials_ts)<330)=[];
[~,idx]=findpeaks(-unit_rot);
idx(unit_rot(idx)>30)=[];

ts=zeros(1,length(trials_ts));
count=1;
for i=1:length(trials_ts)
    count=arrayfun(@(x) trials_ts(i)-x,idx);
    if any(count<2 & count>0)
        ts(i)=1;
    end
end
trials_ts(~ts)=[];
% trials_ts=trials_ts.loc;
% trials_ts(1)=[];