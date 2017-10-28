%% Import positional data from tablet tracker
clear frame_ts

fn='VR20170814182104.csv';
importFlag={'"ca1012"'};
% 
% fn='VR20170813130533.csv';
% importFlag={'"ca1011 start"','"stop 1"'};

fn='VR20170814182104.csv';
importFlag={'"ca1012"','"end"'};


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

idx=[sign(diff(unit_posX));sign(diff(unit_posY))];
distance=zeros(1,size(idx,2));
distance(idx(1,:)==-1 & idx(2,:)==0 | idx(1,:)==1 & idx(2,:)==-1 | idx(1,:)==1 & idx(2,:)==1)=-1;
distance(distance==0)=1;

unit_pos=distance.*sqrt(diff(unit_posX).^2+diff(unit_posY).^2);
cum_pos=cumsum(unit_pos);
cum_pos=[0 cum_pos+cum_pos(1)];
unit_pos=cum_pos;

idx=find(diff(unit_rot)==0);
trials_ts=idx(find(diff(idx)>1))+2;

for i=1:3:length(trials_ts)
    try
        unit_pos(trials_ts(i):trials_ts(i+3)-1)=unit_pos(trials_ts(i):trials_ts(i+3)-1)-unit_pos(trials_ts(i));
    catch
        unit_pos(trials_ts(i):end)=unit_pos(trials_ts(i):end)-unit_pos(trials_ts(i));
    end
end
arms=trials_ts(1:3:length(trials_ts));
unit_pos(1:trials_ts(1)-1)=-unit_pos(1:trials_ts(1)-1)+unit_pos(arms(2)-1);

unit_vel=diff(cum_pos);
unit_vel=[unit_vel(1) unit_vel];
% unit_vel(unit_vel<-10)=0;

