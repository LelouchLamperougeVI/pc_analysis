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

fid=fopen('rsc036.txt');
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
for mouse=1:3
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

for mouse=1:3
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
end
%% Cumulative dists
figure; hold on;
for i=1:length(n)
    [c,centres]=hist(sparsity{i},floor(length(sparsity{i})/4));
    c=cumsum(c./sum(c));

    plot(centres,c,'displayname',windows_list{i});
    xlabel('sparsity (%)');
    ylabel('cumulative freq.');
    axis square
end
ylim([0 1])
legend(gca,'show','location','southeast');

figure; hold on;
for i=1:length(n)
    [c,centres]=hist(SI{i},floor(length(SI{i})/4));
    c=cumsum(c./sum(c));

    plot(centres,c,'displayname',windows_list{i});
    xlabel('spatial info. (bits)');
    ylabel('cumulative freq.');
    axis square
end
ylim([0 1])
legend(gca,'show','location','southeast');

figure; hold on;
for i=1:length(n)
    [c,centres]=hist(width{i},floor(length(width{i})/4));
    c=cumsum(c./sum(c));

    plot(centres,c,'displayname',windows_list{i});
    xlabel('place-field width (cm)');
    ylabel('cumulative freq.');
    axis square
end
ylim([0 1])
legend(gca,'show','location','southeast');


%% High level stats
figure;
[p,tbl,stats] = kruskalwallis(cell2mat(SI),SI_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
% xticklabels({'2','3','4','5','6','7'})
xticklabels(windows_list2)
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
xticklabels(windows_list2)
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
xticklabels(windows_list2)
xlabel('window');
ylabel('PC width mean ranks');

figure;
[p,tbl,stats] = anova1(cell2mat(proportion),proportion_idx,'off');
[c,m]=multcompare(stats,'display','off');
pval=c(:,6);
groups=c(pval<0.05,1:2);
groups=mat2cell(groups,ones(1,size(groups,1)),2);
super_bar_graphs(m(:,1),m(:,2),groups,pval(pval<0.05));
xticklabels({'2','3','4','5','6','7'})
xticklabels(windows_list2)
xlabel('window');
ylabel('Proportion place cells');


%% Window list and other stuff
windows_list={'posterior RCS','anterior RCS','anterior PPC','M2','posterior PPC','S1','M1/S1'};
windows_list2={'pRCS','aRCS','aPPC','M2','pPPC','S1','M1/S1'};

vr_length=150;
bins=50;


%% Sequences across windows

for session=23:36
    for window=1:5
        try
            if ~isempty(analysis(1).session(session).analysis(window).stack)
                    figure;
                    imagesc(analysis(1).session(session).analysis(window).stack(analysis(1).session(session).analysis(window).pc_list,:));
                    colormap jet
                    title([windows_list(str2double(analysis(1).session(session).analysis(window).window(end))) ' ' num2str([session window])]);
                    ylabel('ordered neuron no.');
                    xlabel('position (cm)');
                    set(gca,'xtick',0:bins/5:bins);
                    set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
            end
        catch
        end
    end
end

%% Q matrix
session=33;
window=5;

qMatrix=corr(analysis(1).session(session).analysis(window).stack(analysis(1).session(session).analysis(window).pc_list,:));

figure;
imagesc(qMatrix);
set(gca,'xtick',0:bins/5:bins);
set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
xlabel('position (cm)');
set(gca,'ytick',0:bins/5:bins);
set(gca,'yticklabel',strsplit(num2str(-vr_length:vr_length/5:0)));
ylabel('position (cm)');
c=colorbar; c.Label.String='corr. coef.';
colormap jet
title(windows_list(str2double(analysis(1).session(session).analysis(window).window(end))));
axis square

%% Single place cells
session=33;
window=5;
n=[128 90 112];
% n=[];

list=analysis(1).session(session).analysis(window).pc_list;

count=1;
if isempty(n)
    figure;
    for k=list
        if count>25
            count=1;
            figure;
        end
        subplot(5,5,count);
        imagesc(analysis(1).session(session).analysis(window).psth{k});
        set(gca,'xtick',0:bins/4:bins);
        set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/4:0)));
        title(['n = ' num2str(k)]);
        colormap hot
        ylabel('trials')
        xlabel('distance (cm)')
        colorbar
        count=count+1;
    end
else
    for k=n
        figure
        imagesc(analysis(1).session(session).analysis(window).psth{k});
        set(gca,'xtick',0:bins/4:bins);
        set(gca,'xticklabel',strsplit(num2str(-vr_length:vr_length/4:0)));
        colormap hot
        ylabel('trials')
        xlabel('distance (cm)')
        colorbar
    end
end










