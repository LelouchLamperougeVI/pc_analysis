trunk=dir;
idx=find(vertcat(trunk.isdir));
idx(1:2)=[];
trunk=arrayfun(@(x) [trunk(x).folder '\' trunk(x).name],idx,'uniformoutput',false);

session=[];
count=0;
for i=1:length(trunk)
    branch=dir(trunk{i});
    idx=find(vertcat(branch.isdir));
    idx(1:2)=[];
    branch=arrayfun(@(x) branch(x).name,idx)';
    branch=cellfun(@(x) strfind(branch,num2str(x)), num2cell(1:20),'uniformoutput',false);
    branch=~cellfun(@(x) isempty(x),branch);
    if branch(1:3)
        count=count+1;
        date{count}=trunk{i}(end-9:end);
        path{count}=trunk{i};
        session=[session 1];
    end
    if branch(4:6)
        count=count+1;
        date{count}=trunk{i}(end-9:end);
        path{count}=trunk{i};
        session=[session 2];
    end
    if branch(7:9)
        count=count+1;
        date{count}=trunk{i}(end-9:end);
        path{count}=trunk{i};
        session=[session 3];
    end
end

date=datetime(date,'inputformat','yyyy_MM_dd','format','dd MMMM yyyy');

idx=num2cell(1:9);
ev=zeros(1,length(path));
rev=zeros(1,length(path));
for i=1:length(path)
    switch session(i)
        case 1
            for j=1:3
                load([path{i} '\' num2str(idx{j}) '\behavior.mat']);
                load([path{i} '\' num2str(idx{j}) '\Plane1\deconv.mat']);
                b{j}=behavior;
                d{j}=deconv;
            end
            [ev(i),rev(i)]=explained_variance(d,b,3);
        case 2
            for j=4:6
                load([path{i} '\' num2str(idx{j}) '\behavior.mat']);
                load([path{i} '\' num2str(idx{j}) '\Plane1\deconv.mat']);
                b{j-3}=behavior;
                d{j-3}=deconv;
            end
            [ev(i),rev(i)]=explained_variance(d,b,3);
    end
end

%%
fid=fopen('rsc37_log.csv');
C=textscan(fid,'%s','delimiter',',');
fclose(fid);
C=C{1,1};
C=reshape(C,[],3);
index.date=cellfun(@(x) strrep(x,'"',''),C(:,1),'uniformoutput',false);
index.date=datetime(index.date, 'inputformat','dd/MM/yyyy','format','dd MMMM yyyy');
index.window=cellfun(@(x) str2double(strrep(strrep(x,'"',''),'win','')),C(:,2));
index.session=cellfun(@(x) str2double(strrep(x,'"','')),C(:,3));

window=zeros(1,length(date));
for i=1:length(date)
    window(i)=index.window(logical((date(i)==index.date) .* (session(i)==index.session)));
end

%% Stats
bar_stats([ev' rev'],'groups',window,'g_labels',cellstr(strrep(num2str(0:7),' ','')'),'v_labels',{'ev','rev'},'test','wilcoxon');



bar_stats([ev' rev'],'groups',window,'g_labels',{'pRSC','mRSC','aRSC','M2','M2','PPC','S1','M1'},'v_labels',{'ev','rev'},'test','wilcoxon');








