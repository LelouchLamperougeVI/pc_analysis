function behavior=categories_import
% Import behaviour from object categorization csv

[fn,path]=uigetfile('*.csv','Select CSV file');
cd(path);

fid=fopen(fn);
C=textscan(fid,'%s','delimiter',',');
fclose(fid);
C=C{1,1};
charC=char(C);

c_index=csv_parse2(C);
index=find(c_index==1);
recording_length=zeros(1,length(index));
frame_ts=cell(1,length(index));
for i=1:length(index)-1
    idx=c_index == 2;
    idx([1:index(i)-1 index(i+1)+1:end])=false;
    if strcmp(C(find(idx,1)+4),'0')
        idx(find(idx,1))=false;
    end
    idx=find(idx)-1;
    frame_ts{i}=charArray2double(charC(idx(1:1:end),:)');
    recording_length(i)=length(frame_ts{i});
    fprintf('%d) %d frames\n',i,recording_length(i));
end
idx=c_index == 2;
idx(1:index(end)-1)=false;
if strcmp(C(find(idx,1)+4),'0')
    idx(find(idx,1))=false;
end
idx=find(idx)-1;
frame_ts{end}=charArray2double(charC(idx(1:1:end),:)');
recording_length(end)=length(frame_ts{end});
fprintf('%d) %d frames\n',length(index),recording_length(end));

choice=input('Which recording session to extract?\n');
frame_ts=frame_ts{choice};

choice=input('How many frames to extract? [enter] to extract entire length\n');
if choice
    frame_ts=frame_ts(1:choice);
end

object_ts=[];
black_ts=[];
while(isempty(object_ts) || isempty(black_ts)) %because matlab is stupid...
    idx=find(c_index==4);
    object_ts=charArray_cmp(charC(idx+1,:)','bundle');
    object_type=charC(idx(logical(object_ts))+2,:);
    object_ts=charArray2double(charC(idx(logical(object_ts))-3,:)');
    start=object_ts<frame_ts(1);
    ending=object_ts>frame_ts(end);
    object_ts(start|ending)=[];
    object_ts=frame_ts(knnsearch(frame_ts',object_ts'));
    object_frame=zeros(1,length(frame_ts));
    object_frame(knnsearch(frame_ts',object_ts'))=1;
    object_type(start|ending,:)=[];

    idx=find(c_index==3);
    black_ts=charArray_cmp(charC(idx+2,:)','black');
    black_ts=charArray2double(charC(idx(logical(black_ts))-1,:)');
    start=black_ts<frame_ts(1);
    ending=black_ts>frame_ts(end);
    black_ts(start|ending)=[];
    black_ts=frame_ts(knnsearch(frame_ts',black_ts'));
    black_frame=zeros(1,length(frame_ts));
    black_frame(knnsearch(frame_ts',black_ts'))=1;
end


behavior.frame_ts=frame_ts;
behavior.object_ts=object_ts;
behavior.object_frame=object_frame;
behavior.object_type=object_type;
behavior.black_ts=black_ts;
behavior.black_frame=black_frame;
