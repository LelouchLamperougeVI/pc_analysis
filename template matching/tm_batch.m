dates=dir;
idx=arrayfun(@(x) dates(x).isdir, 1:length(dates));
idx(1:2)=0;
dates=arrayfun(@(x) dates(x).name, find(idx),'uniformoutput',false);

root=pwd;

for i=1:length(dates)
    cd(dates{i});
    
    try
        if exist('2','dir')
            cd 2;
            load behavior.mat;
            load Plane1\deconv.mat;
            load Plane1\timecourses.mat;
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis=pc_batch_analysis(behavior,deconv,'test','mixed','bins',80);

            cd ..\1
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis','analysis','deconv');

            cd ..\3
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis','analysis','deconv');

            cd ..
        end
    catch
        cd(root);
        cd(dates{i});
    end
    
    try
        if exist('5','dir')
            cd 5;
            load behavior.mat;
            load Plane1\deconv.mat;
            load Plane1\timecourses.mat;
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis=pc_batch_analysis(behavior,deconv,'test','mixed','bins',80);

            cd ..\4
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis','analysis','deconv');

            cd ..\6
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis','analysis','deconv');

            cd ..
        end
    catch
        cd(root);
        cd(dates{i});
    end
    
    cd ..
end

%% Reversed replay
dates=dir;
idx=arrayfun(@(x) dates(x).isdir, 1:length(dates));
idx(1:2)=0;
dates=arrayfun(@(x) dates(x).name, find(idx),'uniformoutput',false);

root=pwd;

for i=1:length(dates)
    cd(dates{i});
    
    try
        if exist('2','dir')
            cd 2;
            load behavior.mat;
            load Plane1\deconv.mat;
            load Plane1\timecourses.mat;
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis=pc_batch_analysis(behavior,deconv,'test','mixed','bins',80);

            cd ..\1
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis_reverse','analysis','deconv');

            cd ..\3
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis_reverse','analysis','deconv');

            cd ..
        end
    catch
        cd(root);
        cd(dates{i});
    end
    
    try
        if exist('5','dir')
            cd 5;
            load behavior.mat;
            load Plane1\deconv.mat;
            load Plane1\timecourses.mat;
            [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
            analysis=pc_batch_analysis(behavior,deconv,'test','mixed','bins',80);

            cd ..\4
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis_reverse','analysis','deconv');

            cd ..\6
            load behavior.mat;
            load Plane1\deconv.mat;
            [analysis,deconv]=tm_analysis(analysis,behavior,deconv);
            save('tm_analysis_reverse','analysis','deconv');

            cd ..
        end
    catch
        cd(root);
        cd(dates{i});
    end
    
    cd ..
end