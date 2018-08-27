function behavior=import_sutter(vr_length)

[fn,path]=uigetfile('*.txt');
cd(path);

fid=fopen(fn);
% C=fscanf(fid,'%s,%d,%s,%d,%s,%d\n');
C=textscan(fid,'%s');
fclose(fid);

C=C{1,1};

% C=cellfun(@(x) sscanf(x,'%2c,%f,%2c,%f,%2c,%f'), C,'uniformoutput',false);
C=cellfun(@(x) sscanf(x,'%[tfmlpr,]%f'), C,'uniformoutput',false);
idx=cellfun(@length, C);

temp=C(idx==14);
temp=cell2mat(temp')';

frame_ts=temp(:,4)./1000;
unit_pos=temp(:,14);

temp=C(idx==11);
temp=cell2mat(temp')';

trials=temp(:,4)./1000;
trials_ts=knnsearch(frame_ts,trials);
trials=frame_ts(trials_ts);
trials(trials_ts==1)=[];
trials_ts(trials_ts==1)=[];

for i=trials_ts'
    unit_pos(i:end)=unit_pos(i:end)-unit_pos(i);
end

mean_step=mean(unit_pos(trials_ts(2:end)-1));
max_step=max(unit_pos(trials_ts(2:end)-1));

unit_pos(1:trials_ts(1)-1)=unit_pos(1:trials_ts(1)-1)-unit_pos(trials_ts(1)-1)+mean_step;

unit_pos(1:trials_ts(1)-1)=unit_pos(1:trials_ts(1)-1)./max_step.*vr_length;
for i=1:length(trials_ts)-1
    unit_pos(trials_ts(i):trials_ts(i+1)-1)=unit_pos(trials_ts(i):trials_ts(i+1)-1)./max_step.*vr_length;
end
unit_pos(trials_ts(end):end)=unit_pos(trials_ts(end):end)./max_step.*vr_length;

unit_vel=diff(unit_pos);
unit_vel=[unit_vel(1); unit_vel];
unit_vel(trials_ts)=0;

behavior.frame_ts=frame_ts';
behavior.trials=trials';
behavior.trials_ts=trials_ts';
behavior.unit_pos=unit_pos';
behavior.unit_vel=unit_vel';