function deconv=poisson_sim(varargin)
% Generate simulated firing rates modeled as independent Poisson processes
% Inputs:
%   m: number of samples
%   n: number of neurons
%   'assemblies': index vector of neuronal assemblies e.g. [1 1 0 2 2 0 2]
%   'R': mean spike rate [min max] chosen at random from a uniform
%       distribution between the assigned values (default [1 5])
%   'prob': probability of activation e.g. 0.005
%   'activationR': firing rate during reactivation events (default [20 30])
% Outputs:
%   deconv: you should know what this is...

m=varargin{1};
n=varargin{2};
[assemblies,R,prob,activationR]=parse_input(varargin);

R=ceil((diff(R)+1)*rand(n,1)+R(1)-1);
deconv=arrayfun(@(x) poissrnd(x,m,1),R,'uniformoutput',false);
deconv=reshape(cell2mat(deconv),m,n);

bins=ceil(m*prob);

for i=1:length(unique(assemblies))-1
    events=randi(m,1,bins);
    deconv(events,assemblies==i)=randi(diff(activationR)+1,bins,sum(assemblies==i))+activationR(1)-1;
end

function [assemblies,R,prob,activationR]=parse_input(inputs)
assemblies=[];
R=[1 5];
prob=[];
activationR=[20 30];

idx=3;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'assemblies'
            idx=idx+1;
            assemblies=inputs{idx};
        case 'r'
            idx=idx+1;
            R=inputs{idx};
        case 'prob'
            idx=idx+1;
            prob=inputs{idx};
        case 'activationr'
            idx=idx+1;
            prob=inputs{idx};
        otherwise
    end
    idx=idx+1;
end
