function [analysis,deconv]=tm_analysis(analysis,behavior,deconv,cf_range,shuffles)
% Template matching analysis
% Generates template from behavior and deconv included in analysis output
% from pc_batch_analysis
% Matched with pre/post recordings with cell identities shuffling test
% Inputs:
%   'analysis': from pc_batch_analysis
%   'behavior': pre/post-task behavior
%   'deconv':   pre/post-task deconv
%   'cf_range': range of compression factors to test (default - [1 20])
%               can also be vector with values of desired CF's to be tested
%   'shuffles': number of identity shuffles (default - 1,000)

if nargin<4
    cf_range=[1 20];
end
if nargin<5
    shuffles=100;
end
if length(cf_range)==2
    cf_range=cf_range(1):cf_range(2);
elseif length(cf_range)<2
    error('cf_range must be a vector of positive integers with at least 2 values');
end

template=get_template(analysis);
deconv=noRun_deconv(behavior,deconv);

%reversed replay
% template=template(end:-1:1,:);

count=0;
f=waitbar(count,'Template matching...');
C=NaN(size(deconv,1),length(cf_range));
pval=C;
for i=cf_range(end:-1:1)
    t=compress(template,i);
    [temp,p]=template_match(t,deconv,shuffles);
    idx=floor((size(C,1)-length(temp))/2);
    C(idx+1:idx+length(temp),i)=temp;
    pval(idx+1:idx+length(temp),i)=p;
    count=count+1;
    waitbar(count/length(cf_range),f);
end
close(f);

analysis.template=t;
analysis.C=C;
analysis.C_pval=pval;



function t=compress(template,i)
% this method should be faster???
kernel=(1/i).*ones(i,size(template,2)); %averaging kernel
t=mat_conv(template,kernel);
idx=ceil(i/2):i:size(t,1);
t=t(idx,:);



function deconv=noRun_deconv(behavior,deconv)
try
    unit_vel=behavior.unit_vel;
catch
    unit_vel=behavior.speed_raw;
end
try
    thres=noRun(unit_vel);
    thres=abs(unit_vel)<thres;
    deconv=deconv(thres,:);
catch
    return;
end