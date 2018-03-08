function behavior=tablet_import
% Import positional data from tablet tracker

[fn,path]=uigetfile('*');
cd(path);

fid=fopen(fn);
C=textscan(fid,'%s','delimiter',',');
fclose(fid);
C=C{1,1};

c_index=csv_parse(C);
index=find(c_index==1);
idx=list_box(['start';C(index+1);'end']);
if idx{1}==1
    recStart=1;
else
    recStart=index(idx{1}-1);
end
if idx{2}==length(C(index+1))+2
    recFinish=length(C);
else
    recFinish=index(idx{2}-1);
end

if recStart>=recFinish
    error('Stop messing around you dummy!')
end

c_index(1:recStart)=0;
c_index(recFinish:end)=0;

C=char(C);

pos=C(find(c_index==2)+2,:);
pos=charArray2double(pos');

ts=C(find(c_index==2)-1,:);
ts=charArray2double(ts');

frame_ts=C(find(c_index==3)-1,:);
flag=C(find(c_index==3)+2,:);
if strcmp(flag(1,1),'1')
    frame_ts(1,:)=[];
end
if strcmp(flag(end,1),'1')
    frame_ts(end)=[];
end
frame_ts=charArray2double(frame_ts');
frame_ts=frame_ts(1:2:end);

trials=C(find(c_index==4)+1,:);
idx=charArray_cmp(trials','trial');
trials=C(find(c_index==4)+2,:);
idx=logical(idx) & logical(charArray_cmp(trials','trigger'));
trials=C(find(c_index==4)-1,:);
trials=charArray2double(trials(idx,:)');

trials(trials>frame_ts(end) | trials<frame_ts(1))=[];

%make unit_pos, cum_pos and unit_vel
try
    unit_pos=arrayfun(@(x) find(ts>=x,1), frame_ts);
    unit_pos=pos(unit_pos);
catch
    ts=[ts frame_ts(end)];
    pos=[pos pos(end)];
    unit_pos=arrayfun(@(x) find(ts>=x,1), frame_ts);
    unit_pos=pos(unit_pos);
end

idx=find(diff(unit_pos)<-10);
cum_pos=unit_pos(1:idx(1));
for i=1:length(idx)-1
    cum_pos=[cum_pos unit_pos(idx(i)+1:idx(i+1))-unit_pos(idx(i)+1)+cum_pos(end)];
end
cum_pos=[cum_pos unit_pos(idx(end)+1:end)-unit_pos(idx(end)+1)+cum_pos(end)];

unit_vel=diff(cum_pos)./diff(frame_ts);
unit_vel=[unit_vel(1) unit_vel];
unit_vel(isnan(unit_vel))=0;
% unit_vel(unit_vel<-10)=0;

% sync trials index
trials_ts=arrayfun(@(x) find(frame_ts>=x,1),trials);

% trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)=trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)+1;
% trials=frame_ts(trials_ts);

% trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)=trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)+1;
% trials=frame_ts(trials_ts);

trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)=trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)+1;
trials=frame_ts(trials_ts);

% compatibility with Dun's code
behavior.ts=ts;
behavior.pos_raw=pos;
behavior.pos_norm=(pos-min(pos))./max(pos-min(pos));
idx=find(diff(pos)<-10);
behavior.pos_cum=pos(1:idx(1));
for i=1:length(idx)-1
    behavior.pos_cum=[behavior.pos_cum pos(idx(i)+1:idx(i+1))-pos(idx(i)+1)+behavior.pos_cum(end)];
end
behavior.pos_cum=[behavior.pos_cum pos(idx(end)+1:end)-pos(idx(end)+1)+behavior.pos_cum(end)];
behavior.speed=diff(behavior.pos_cum)./diff(ts);
behavior.speed=[behavior.speed(1) behavior.speed];
idx=isnan(behavior.speed) | isinf(behavior.speed);
behavior.speed(idx)=[];
ts(idx)=[];
behavior.speed=interp1(ts,behavior.speed,behavior.ts);
behavior.speed_raw=unit_vel';
behavior.trial=bsxfun(@ge,trials',ts);
behavior.trial=sum(behavior.trial);
behavior.trial=abs(behavior.trial-max(behavior.trial))+1;
% for my code
if mod(length(frame_ts),10)
    error('Uneven frame_ts. Direct complaint to HaoRan');
end
behavior.frame_ts=frame_ts;
behavior.trials=trials;
behavior.trials_ts=trials_ts;
behavior.unit_pos=unit_pos;
behavior.unit_vel=unit_vel;


