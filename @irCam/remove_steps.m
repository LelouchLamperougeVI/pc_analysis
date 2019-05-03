function remove_steps(obj,thres,samples)
%detect abrupt step changes in signal baseline and remove them
%thres: fraction of total residual for minimum stepping (default 0.05)
%samples: minimum number of samples a step contains (default 10 sec)

if isempty(obj.traces)
    error('must extract traces first');
end

if nargin<2
    thres=.1; % 10% of the max decrease in residual
end
if nargin<3
    samples = obj.cam.FrameRate * 10; %10 second gaps
end

baseline = obj.original_traces(:,obj.get_label('baseline'));

residue0 = length(baseline) * log(var(baseline));
[~,residue1] = cpsingle(baseline, 'std', samples);

thres = (residue0 - residue1) * thres;

steps=findchangepts(baseline,'statistic','std', 'minthreshold',thres, 'mindistance',samples);
steps=[1; steps];

temp=[];
for ii=1:length(steps)-1
    temp = [temp; zscore(diff(obj.original_traces(steps(ii):steps(ii+1)-1, :)))];
%     obj.traces(steps(ii):steps(ii+1)-1, :)=obj.original_traces(steps(ii):steps(ii+1)-1, :) - mean( obj.original_traces(steps(ii):steps(ii+1)-1, :) );
end
% obj.traces(steps(end):end, :)=obj.original_traces(steps(end):end, :) - mean( obj.original_traces(steps(end):end, :) );
temp = [temp; zscore(diff(obj.original_traces(steps(end):end, :)))];

obj.traces = temp;