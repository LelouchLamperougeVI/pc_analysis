function behavior=import_2d(fn,path)
%import behavior from old 2d tablet navigation csv files

if ~exist('fn','var')
    [fn,path]=uigetfile('*.csv');
end
cd(path);

fid=fopen(fn);
C=textscan(fid,'%s');
fclose(fid);

C=C{1,1};

x=[];y=[];ts=[];reward=[];punish=[];trial=[];objects=[];obj_probe=[];
for i=1:length(C)
    temp=textscan(C{i},'%.4f,position,%.4f,%.4f');
    if ~isempty(temp{2})
        x=[x;temp{2}];
        y=[y;temp{3}];
        ts=[ts;temp{1}];
    end
    temp=textscan(C{i},'%.4f,reward,%.4f');
    if ~isempty(temp{2})
        reward=[reward;temp{1}];
    end
    temp=textscan(C{i},'%.4f,airPuff,%.4f');
    if ~isempty(temp{2})
        punish=[punish;temp{1}];
    end
    temp=textscan(C{i},'%.4f,trial,%.4f');
    if ~isempty(temp{2})
        trial=[trial;temp{1}];
    end
    
    temp=textscan(C{i},'%.4f,%[objects,cone]%s');
    if isempty(temp{3}) && ~isempty(temp{2})
        objects=[objects i];
        obj_probe=[obj_probe; temp{1}];
    end
end
trial=knnsearch(ts,trial);
trial=unique(trial);
reward=knnsearch(ts,reward);
punish=knnsearch(ts,punish);

if length(objects)>1
    behavior.probe=true;
else
    behavior.probe=false;
end

cone=cell(1,length(objects));
sphere=cone;
for i=1:length(objects)
    count=objects(i);
    cone{i}=[str2double(C{count+1}) str2double(C{count+2})];
    sphere{i}=[];
    count=count+12;
    while strcmp(C{count},'cone,') || strcmp(C{count},'sphere,')
        if strcmp(C{count},'cone,')
            cone{i}=[cone{i}; str2double(C{count+1}) str2double(C{count+2})];
        elseif strcmp(C{count},'sphere,')
            sphere{i}=[sphere{i}; str2double(C{count+1}) str2double(C{count+2})];
        end
        count=count+12;
    end
end

if length(objects)==1
    cone=cone{1,1};
    sphere=sphere{1,1};
else
    count=1;
    while count<length(obj_probe)
        if sum(ts(trial)>obj_probe(count) & ts(trial)<obj_probe(count+1))<2
            obj_probe(count)=[];
            cone(count)=[];
            sphere(count)=[];
        else
            count=count+1;
        end
    end
end

temp=split(path,'\');
behavior.mouse=temp{end-1};

temp=textscan(fn,'VR%8d');
behavior.date=temp{1,1};

behavior.x=x;
behavior.y=y;
behavior.ts=ts;
behavior.reward=reward;
behavior.punish=punish;
behavior.trial=trial;
behavior.cone=cone;
behavior.sphere=sphere;
behavior.obj_probe=obj_probe;