function remove_steps(obj,thres,samples)
%detect abrupt step changes in signal baseline and remove them
%thres: fraction of total residual for minimum stepping (default 0.05)
%samples: minimum number of samples a step contains (default 10 sec)

if isempty(obj.traces)
    error('must extract traces first');
end

if nargin<2
    thres=.05;
end
if nargin<3
    samples = 10;
end

samples=round(samples * obj.cam.FrameRate);
base_residuals=sum(abs(obj.traces-mean(obj.traces)));

steps=findchangepts(obj.traces,'statistic','mean','minthreshold',base_residuals*thres,'mindistance',samples);
steps=[1 steps];

for ii=1:length(steps)-1
    obj.traces(steps(ii):steps(ii+1)-1)=obj.traces(steps(ii):steps(ii+1)-1)-mean(obj.traces(steps(ii):steps(ii+1)-1));
end
obj.traces(steps(end):end)=obj.traces(steps(end):end)-mean(obj.traces(steps(end):end));