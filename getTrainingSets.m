count=1;
sigh=dir;
for l=3:length(sigh)
    progressbar((l-3)/(length(sigh)-3));
    cd(sigh(l).name);
    d=dir;
    for i=3:length(d)
        cd(d(i).name);
        dd=dir;
        for j=3:length(dd)
            try
                clear dat
                cd(dd(j).name);
                f=dir('*_proc.mat');
                load(f.name);
                f=dir('*_Nk300.mat');
                load(f.name);
%                 [Input,Mimg,Target,features]=genSets(dat,Fcell);
%                 [Target,features]=genSetsFeatures(dat,Fcell);
%                 [Input,Target]=genSetsFeatures(dat,Fcell);
                [Input,Target]=genSetsInput(dat,Fcell);
                I{count}=Input;
%                 M{count}=Mimg;
                T{count}=Target;
%                 F{count}=features;
                count=count+1;
            catch
            end
            cd ..
        end
        cd ..
    end
    cd ..
end

%%
Input=[];
Mimg=[];
Target=[];
Features=[];
for i=1:length(T)
    Input=[Input I{i}];
%     Mimg=[Mimg M{i}];
    Target=[Target T{i}];
%     Features=[Features F{i}];
end