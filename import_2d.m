function behavior=import_2d(plotFlag)
%import behavior from old 2d tablet navigation csv files

if nargin<1
    plotFlag=false;
end

[fn,path]=uigetfile('*.csv');
cd(path);

fid=fopen(fn);
% C=fscanf(fid,'%s,%d,%s,%d,%s,%d\n');
C=textscan(fid,'%s');
fclose(fid);

C=C{1,1};

x=[];y=[];ts=[];reward=[];punish=[];trial=[];
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
end
trial=knnsearch(ts,trial);
reward=knnsearch(ts,reward);
punish=knnsearch(ts,punish);

found=false;count=1;
while ~found
    temp=textscan(C{count},'%.4f,%[objects,cone]%s');
    if isempty(temp{3}) && ~isempty(temp{2})
        found=true;
    else
        count=count+1;
    end
end
cone=[str2double(C{count+1}) str2double(C{count+2})];
sphere=[];
count=count+12;
while strcmp(C{count},'cone,') || strcmp(C{count},'sphere,')
    if strcmp(C{count},'cone,')
        cone=[cone; str2double(C{count+1}) str2double(C{count+2})];
    elseif strcmp(C{count},'sphere,')
        sphere=[sphere; str2double(C{count+1}) str2double(C{count+2})];
    end
    count=count+12;
end

if plotFlag
    figure;
    hold on;
    plot(cone(:,1),cone(:,2),'g^','markerfacecolor','g','markersize',20);
    plot(sphere(:,1),sphere(:,2),'ro','markerfacecolor','r','markersize',20);
    
    xlim([-100 100]);
    ylim([0 140]);
    
    for i=randperm(length(trial),10) %randomly plot 10 trials
        plot(x(trial(i)+1:trial(i+1)),y(trial(i)+1:trial(i+1)));
        idx=reward(reward>trial(i) & reward<=trial(i+1));
        plot(x(idx),y(idx),'kd','markersize',20);
        idx=punish(punish>trial(i) & punish<=trial(i+1));
        plot(x(idx),y(idx),'bx','markerfacecolor','b','markersize',20);
    end
end

temp=split(path,'\');
behavior.mouse=temp{end-1};

behavior.date=textscan(fn,'VR%8d');

behavior.x=x;
behavior.y=y;
behavior.ts=ts;
behavior.reward=reward;
behavior.punish=punish;
behavior.trial=trial;
behavior.cone=cone;
behavior.sphere=sphere;