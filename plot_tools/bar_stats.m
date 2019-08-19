function h = bar_stats(varargin)
% Plot stats (still need to implement multiple comparisons)
%
% Inputs
%   x:  matrix with R observations and C variables
%   'groups' (optional): vector index of groups
%   'g_labels': group labels
%   'v_labels': variable labels
%   'type':
%       'mean' (default for ttest)
%       'median' (default for wilcoxon)
%   'test' (optional):
%       'ttest'
%       'wilcoxon'

x=varargin{1};
[groups,type,testFlag,g_labels,v_labels]=parse_input(varargin);
p=zeros(1,length(unique(groups)));
count=1;
for i=unique(groups)
    idx=groups==i;
    X=zeros(sum(idx),size(x,2));
    for k=1:size(x,2)
        X(:,k)=x(idx,k);
    end
    if testFlag==1
        [~,p(count)]=ttest2(X(:,1),X(:,2));
    elseif testFlag==2
        [p(count),~,stats]=ranksum(X(:,1),X(:,2));
    end
    if type==1
        line=mean(X);
    elseif type==2
        line=median(X);
    end
    count=count+1;
end
p(p>0.05)=NaN;
h=grouped_boxplot(x,groups,g_labels,v_labels,p);


function [groups,type,testFlag,g_labels,v_labels]=parse_input(inputs)
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
        case 'g_labels'
            idx=idx+1;
            g_labels=inputs{idx};
        case 'v_labels'
            idx=idx+1;
            v_labels=inputs{idx};
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