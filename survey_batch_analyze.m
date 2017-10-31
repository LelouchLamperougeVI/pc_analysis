%% Batch analysis for place cells survey

%% Get recording structures
folder='X:\scratch\';
folder='Y:\homes\dun.mao\2photon\raw\';
folder='X:\homes\ingrid.esteves\2P\';
mice={'rsc036','rsc037','rsc038','rsc039'};

cd(folder);

for i=1:length(mice)
    survey(i).mouse=mice(i);
    
    cd(mice{i});
    l1=dir;
    for j=3:length(l1)
        try
            cd(l1(j).name);

            try
                cd 2
                notes=xml2struct('Experiment.xml');
                notes=notes.ThorImageExperiment.ExperimentNotes.Attributes.text;
                notes=split(notes);
                notes=notes{1};
                survey(i).session(j-2).window{1}=notes;
                survey(i).session(j-2).date=l1(j).name;
                cd ..
            catch
            end

            try
                cd 5
                notes=xml2struct('Experiment.xml');
                notes=notes.ThorImageExperiment.ExperimentNotes.Attributes.text;
                notes=split(notes);
                notes=notes{1};
                survey(i).session(j-2).window{2}=notes;
                survey(i).session(j-2).date=l1(j).name;
                cd ..
            catch
            end

            cd ..
        catch
        end
    end
    
    cd ..
end

%% Batch analyze data

folder='X:\homes\ingrid.esteves\analysis\';
mice={'rsc036','rsc037','rsc038'};

cd(folder);

for i=1:length(mice)
    disp(['Starting mouse #' num2str(i)])
    
    analysis(i).mouse=mice(i);
    
    cd(mice{i});
    l1=dir;
    for j=3:length(l1)
        disp([mice{i} ' ' num2str(j-2) '/' num2str(length(l1)-2)])
        
        analysis(i).session(j-2).date=l1(j).name;
        
        cd(l1(j).name);
        l2=dir;
        
        for k=3:length(l2)
            if l2(k).isdir
                cd(l2(k).name);
                try
                    load behavior.mat
                    load Plane1\deconv.mat
                    load Plane1\timecourses.mat
                    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
                    analysis(i).session(j-2).analysis(k-2)=pc_batch_analysis(behavior,deconv);
                    analysis(i).session(j-2).nRecordings=str2double(l2(k).name);
                catch
                end
                cd ..
            end
        end
        cd ..
    end
    cd ..
end



%% Get index structure using image logs

fid=fopen('rsc038.txt');
C=textscan(fid,'%s','delimiter',',');
fclose(fid);
C=C{1,1};

idx=cellfun(@(x) strcmp(x, {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'}), C,'uniformoutput',false);
idx=cell2mat(idx);

idx1=find(idx(:,1));
for i=1:sum(idx(:,1))
    index(i).date=C{idx1(i)-1};
    try
        idx2=idx;
        idx2(1:idx1(i)-1,:)=0;
        idx2(idx1(i+1):end,:)=0;
    catch
        idx2=idx;
        idx2(1:idx1(i)-1,:)=0;
    end
    index(i).session=find(any(idx2));
    index(i).param=C(find(any(idx2,2))+1);
    index(i).window=C(find(any(idx2,2))+2);
end

%% Add window tag to analysis
mouse=1;
for i=1:length(analysis(mouse).session)
    date=analysis(mouse).session(i).date;
    idx=find(strcmp({index(mouse).index(:).date},strrep(date,'_','/')));
    if ~isempty(idx)
        for j=1:length(analysis(mouse).session(i).analysis)
            try
                analysis(mouse).session(i).analysis(j).window=index(mouse).index(idx).window{j}(1:4);
                analysis(mouse).session(i).analysis(j).param=index(mouse).index(idx).param{j};
            catch
            end
        end
    end
end

%% Compare PC features across windows
mouse=1;

sparsity=cell(1,7);
SI=cell(1,7);
width=cell(1,7);
proportion=cell(1,7);
n=zeros(1,7);
sparsity_idx=[];
SI_idx=[];
width_idx=[];
proportion_idx=[];

for i=1:length(analysis(mouse).session)
    for j=1:length(analysis(mouse).session(i).analysis)
        try
            if strcmpi(analysis(mouse).session(i).analysis(j).param,'RUN')
                switch analysis(mouse).session(i).analysis(j).window
                    case 'win1'
                        plotID=1;
                    case 'win2'
                        plotID=2;
                    case 'win3'
                        plotID=3;
                    case 'win4'
                        plotID=4;
                    case 'win5'
                        plotID=5;
                    case 'win6'
                        plotID=6;
                    case 'win7'
                        plotID=7;
                    otherwise
                        plotID=0;
                end

                if plotID
                    n(plotID)=n(plotID)+1;
                    sparsity{plotID}=[sparsity{plotID} analysis(mouse).session(i).analysis(j).sparsity];
                    SI{plotID}=[SI{plotID} analysis(mouse).session(i).analysis(j).SI];
                    width{plotID}=[width{plotID} analysis(mouse).session(i).analysis(j).width];
                    proportion{plotID}=[proportion{plotID} sum(analysis(mouse).session(i).analysis(j).pval)/length(analysis(mouse).session(i).analysis(j).pval)];
                    
                    sparsity_idx=[sparsity_idx repmat(plotID,1,length(analysis(mouse).session(i).analysis(j).sparsity))];
                    SI_idx=[SI_idx repmat(plotID,1,length(analysis(mouse).session(i).analysis(j).SI))];
                    width_idx=[width_idx repmat(plotID,1,length(analysis(mouse).session(i).analysis(j).width))];
                    proportion_idx=[proportion_idx plotID];
                end
            end
        catch
        end
    end
end

figure;
[p,tbl,stats] = kruskalwallis(cell2mat(SI),SI_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
xticklabels({'2','3','4','5','6','7'})
xlabel('window');
ylabel('Spatial info mean ranks');

figure;
[p,tbl,stats] = kruskalwallis(cell2mat(sparsity),sparsity_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
xticklabels({'2','3','4','5','6','7'})
xlabel('window');
ylabel('Sparsity mean ranks');

figure;
[p,tbl,stats] = kruskalwallis(cell2mat(width),width_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
xticklabels({'2','3','4','5','6','7'})
xlabel('window');
ylabel('PC width mean ranks');

figure;
[p,tbl,stats] = kruskalwallis(cell2mat(proportion),proportion_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
xticklabels({'2','3','4','5','6','7'})
xlabel('window');
ylabel('Proportion place cells mean ranks');











