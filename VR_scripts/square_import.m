%% Import positional data from tablet tracker
clear frame_ts

fn='VR20170814182104.csv';
importFlag={'"ca1012"'};

fn='VR20170813145900.csv';
importFlag={'"ca1012 start 1"','"stop"'};
% 
% fn='VR20170814185005.csv';
% importFlag={'"ca1011"'};


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

idx=cellfun(@(x) strcmp(x, 'position'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
for i=1:length(idx)
    count=C{idx(i)+1};
    posX(i)=str2double(count);
    count=C{idx(i)+2};
    posY(i)=str2double(count);
    count=C{idx(i)-1};
    ts(i)=str2double(count);
end

idx=cellfun(@(x) strcmp(x, 'rotation'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
for i=1:length(idx)
    count=C{idx(i)+3};
    rot(i)=str2double(count);
    count=C{idx(i)-1};
    rot_ts(i)=str2double(count);
end

idx=cellfun(@(x) strcmp(x, 'pickup'), C);
idx(1:recStart)=0;
if length(importFlag)==2
    idx(recFinish:end)=0;
end
idx=find(idx);
count=1;
for i=1:length(idx)
    if strcmp(C{idx(i)+1},'tunnel') & strcmp(C{idx(i)+2},'enter')
        tunnel(count)=str2double(C{idx(i)-1});
        count=count+1;
    end
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

%make unit_pos, cum_pos and unit_vel
idx=arrayfun(@(x) find(ts>=x,1), frame_ts);
unit_posX=posX(idx);
unit_posY=posY(idx);
unit_posX=unit_posX-min(unit_posX);
unit_posY=unit_posY-min(unit_posY);
idx=arrayfun(@(x) find(rot_ts>=x,1), frame_ts);
unit_rot=rot(idx);

idx=find(diff(unit_rot)==0);
trials_ts=idx(find(diff(idx)>1))+2;

if unit_posX(trials_ts(1))-unit_posX(trials_ts(2))>10
    alt=0;
    unit_pos=-unit_posY(1:trials_ts(1)-1);
elseif unit_posX(trials_ts(1))-unit_posX(trials_ts(2))<10
    alt=2;
    unit_pos=unit_posY(1:trials_ts(1)-1);
elseif unit_posY(trials_ts(1))-unit_posY(trials_ts(2))>10
    alt=1;
    unit_pos=unit_posX(1:trials_ts(1)-1);
else
    alt=3;
    unit_pos=-unit_posX(1:trials_ts(1)-1);
end

for i=1:length(trials_ts)-1
    switch alt
        case 0
            unit_pos=[unit_pos unit_pos(end)-unit_posX(trials_ts(i):trials_ts(i+1)-1)+unit_posX(trials_ts(i))];
        case 1
            unit_pos=[unit_pos unit_pos(end)-unit_posY(trials_ts(i):trials_ts(i+1)-1)+unit_posY(trials_ts(i))];
        case 2
            unit_pos=[unit_pos unit_pos(end)+unit_posX(trials_ts(i):trials_ts(i+1)-1)];
        case 3
            unit_pos=[unit_pos unit_pos(end)+unit_posY(trials_ts(i):trials_ts(i+1)-1)];
    end
    alt=mod(alt+1,4);
end
switch alt
    case 0
        unit_pos=[unit_pos unit_pos(end)-unit_posX(trials_ts(end):end)+unit_posX(trials_ts(end))];
    case 1
        unit_pos=[unit_pos unit_pos(end)-unit_posY(trials_ts(end):end)+unit_posY(trials_ts(end))];
    case 2
        unit_pos=[unit_pos unit_pos(end)+unit_posX(trials_ts(end):end)];
    case 3
        unit_pos=[unit_pos unit_pos(end)+unit_posY(trials_ts(end):end)];
end
cum_pos=unit_pos;
for i=1:4:length(trials_ts)
    try
        unit_pos(trials_ts(i):trials_ts(i+4)-1)=unit_pos(trials_ts(i):trials_ts(i+4)-1)-unit_pos(trials_ts(i));
    catch
        unit_pos(trials_ts(i):end)=unit_pos(trials_ts(i):end)-unit_pos(trials_ts(i));
    end
end
arms=trials_ts(1:4:length(trials_ts));
unit_pos(1:trials_ts(1)-1)=-unit_pos(1:trials_ts(1)-1)+unit_pos(arms(2)-1);

unit_vel=diff(cum_pos);
unit_vel=[unit_vel(1) unit_vel];
% unit_vel(unit_vel<-10)=0;
tunnel(tunnel>frame_ts(end))=[];
tunnel=arrayfun(@(x) find(frame_ts>=x,1),tunnel);

%% Plot thingy
subplot(4,1,[3 4]);
plot(1:length(unit_pos),unit_pos);

%% Export to "behavior" and tcs.tt
behavior.ts=ts;
behavior.pos_raw=pos-min(pos);
behavior.pos_norm=behavior.pos_raw./max(behavior.pos_raw);
behavior.trial=[];
behavior.pos_cum=[];
pre=0;
add_cum=0;
for i=1:length(trials)
    next=find(ts>=trials(i),1);
    behavior.trial=[behavior.trial i.*ones(1,length(ts(pre+1:next)))];
    behavior.pos_cum=[behavior.pos_cum pos(pre+1:next)+add_cum];
    add_cum=pos(next-1);
    pre=next;
end
behavior.trial=[behavior.trial (i+1).*ones(1,length(ts(pre+1:end)))];
behavior.pos_cum=[behavior.pos_cum pos(pre+1:end)+add_cum];
behavior.speed=diff(behavior.pos_cum);
behavior.speed=[behavior.speed behavior.speed(end)];

tcs.tt=frame_ts;


