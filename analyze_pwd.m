function analysis = analyze_pwd(varargin)

load('behavior.mat');
og_behavior=behavior;

fs=struct2cell(dir);
fs=fs(1,:);

numPlanes=0;
for i=1:length(fs)
    numPlanes=numPlanes+double(~isempty(strfind(fs{i},'plane')));
end

for i=1:numPlanes
    cd([pwd '\plane' num2str(i)]);
    load('deconv.mat');
    load('timecourses.mat');
    load('mean_img.mat');
    load('masks_neurons.mat');
    
    behavior=og_behavior;
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    
    analysis(i)=pc_batch_analysis(behavior,deconv,'mask', maskNeurons, mimg,'test','mixed');
    cd ..
end

analysis=merge_planes_fcn(analysis);

function analysis=merge_planes_fcn(analysis)
while length(analysis) > 1
    analysis(end-1)=merge_planes(analysis(end-1),analysis(end));
    analysis(end)=[];
end