%% Import positional data from tablet tracker
clear frame_ts

if ~exist('varset')
    varset=who;
end

fn='VR20170925143214.csv';
importFlag={'"ca1011 tab 2"'};
% 
% fn='VR20170820163035.csv';
% importFlag={'"ca1012 classroom w cues belt"'};

fn='VR20170924160716.csv';
importFlag={'"ca1011 tread"'};

fn='VR20170925143214.csv';
importFlag={'"ca1011 tread"'};

fn='VR20170926185442.csv';
% importFlag={'"ca1011 tread"','"done"'};
importFlag={'"another"'};

% fn='VR20170927153204.csv';
% importFlag={'"go"'};

fn='VR20170930191750.csv';
importFlag={'"ca1011rec2"'};
% fn='VR20170910114202[1].csv';
% importFlag={'"ca1011 tab 2"'};

% fn='VR20171007131001.csv';
% importFlag={'"another"'};

fn='VR20170910114202[1].csv';
importFlag={'"ca1011 tab 2"'};

fn='VR20171025094004.csv';
importFlag={'"rsc36 win1 no tunnel"','"stoped 1"'};




if ~exist('C')
    fid=fopen(fn);
    C=textscan(fid,'%s','delimiter',',');
    fclose(fid);
    C=C{1,1};
end

index=cellfun(@(x) strcmp(x, 'userEntry'), C);
index=find(index);
idx=list_box(['start';C(index+1);'end']);
if idx{1}==1
    recStart=1;
else
    recStart=index(idx{1}+1);
end
if idx{2}==length(C(index+1))+2
    recFinish=length(C);
else
    recFinish=index(idx{2}+1);
end


% recStart=cellfun(@(x) strcmp(x,importFlag{1}), C);
% recStart=find(recStart,1);

% if length(importFlag)==2
%     recFinish=cellfun(@(x) strcmp(x,importFlag{2}), C);
%     recFinish=find(recFinish,1);
% end


index=cellfun(@(x) strcmp(x, {'position','count-21','pickup'}), C,'uniformoutput',false);
index=cell2mat(index);
index(1:recStart,:)=0;
if length(importFlag)==2
    index(recFinish:end,:)=0;
end

idx=find(index(:,1));
for i=1:length(idx)
    count=C{idx(i)+2};
    pos(i)=str2double(count);
    count=C{idx(i)-1};
    ts(i)=str2double(count);
end

idx=find(index(:,2));
for i=1:length(idx)
    count=C{idx(i)-1};
    frame_ts(i)=str2double(count);
end
frame_ts=frame_ts(1:2:end);

idx=find(index(:,3));
for i=1:length(idx)
    count=C{idx(i)-1};
    if strcmp(C{idx(i)+1},'trial')
        if strcmp(C{idx(i)+2},'trigger')
            trials(i)=str2double(count);
        else
            trials(i)=0;
        end
    else
        trials(i)=0;
    end
end
trials(trials==0)=[];
trials(trials>frame_ts(end) | trials<frame_ts(1))=[];

%make unit_pos, cum_pos and unit_vel
unit_pos=arrayfun(@(x) find(ts>=x,1), frame_ts);
unit_pos=pos(unit_pos);

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

trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)=trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)+1;
trials=frame_ts(trials_ts);

trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)=trials_ts(unit_pos(trials_ts)>min(unit_pos)./2)+1;
trials=frame_ts(trials_ts);

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

    
varset=setdiff(who,varset);
varset=setdiff(varset,{'C','frame_ts','trials','trials_ts','cum_pos','unit_pos','unit_vel','behavior','varset'});
clear(varset{:});