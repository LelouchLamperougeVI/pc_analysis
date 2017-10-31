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











