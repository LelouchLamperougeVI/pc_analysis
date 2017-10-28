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
mice={'rsc036','rsc037','rsc038','rsc039'};

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
        
        if exist('2'); cd 2; goback=true; else; goback=false; end
        try
            load behavior.mat
            load Plane1\deconv.mat
            load Plane1\timecourses.mat
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis(i).session(j-2).analysis(1)=pc_batch_analysis(behavior,deconv);
            analysis(i).session(j-2).nRecordings=1;
        catch
        end
        if goback; cd ..; end
        
        if exist('5'); cd 5; goback=true; else; goback=false; end
        try
            load behavior.mat
            load Plane1\deconv.mat
            load Plane1\timecourses.mat
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis(i).session(j-2).analysis(2)=pc_batch_analysis(behavior,deconv);
            analysis(i).session(j-2).nRecordings=2;
        catch
        end
        if goback; cd ..; end
        
        cd ..
    end
    
    cd ..
end




















