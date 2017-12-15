function h = bar_stats(varargin)
% Plot stats (still need to implement multiple comparisons)
%
% Inputs
%   x:  matrix with R observations and C variables
%   'groups' (optional): vector index of groups
%   'type':
%       'mean' (default for ttest)
%       'median' (default for wilcoxon)
%   'test' (optional):
%       'ttest'
%       'wilcoxon'

x=varargin{1};
[groups,type,testFlag]=parse_input(varargin);
for i=unique(groups)
    idx=groups==i;
    X=zeros(sum(idx),size(x,2));
    for k=1:size(x,2)
        X(:,k)=x(idx,k);
    end
    if testFlag==1
        [~,p]=ttest2(X(:,1),X(:,2));
    elseif testFlag==2
        [p,~,stats]=ranksum(X(:,1),X(:,2));
    end
    if type==1
        line=mean(X);
    elseif type==2
        line=median(X);
    end
end


function [groups,type,testFlag]=parse_input(inputs)
groups=1:size(inputs{1});
type=1;
testFlag=0;
typecast=0;

idx=2;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'test'
            idx=idx+1;
            switch inputs{idx}
                case 'ttest'
                    testFlag=1;
                case 'wilcoxon'
                    testFlag=2;
                    if ~typecast
                        type=2;
                    end
                otherwise
                    error('not a valid test');
            end
        case 'groups'
            idx=idx+1;
            groups=inputs{idx};
        case 'type'
            idx=idx+1;
            switch inputs{idx}
                case 'mean'
                    type=1;
                case 'median'
                    type=2;
                otherwise
                    error('not a valid type');
            end
            typecast=1;
        otherwise
    end
    idx=idx+1;
end