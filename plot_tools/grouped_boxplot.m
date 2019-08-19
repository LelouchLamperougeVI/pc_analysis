function h = grouped_boxplot(x,groups,g_labels,v_labels,pval)
% Grouped box plot with multiple variables
% Inputs
%   x: matrix with R observations and C variables
%   groups: grouping numeric vector
%   g_labels: group labels
%   v_labels: variable labels
%   pval: matrix with R between-pairs p-val and C groups

colors_idx='rbkmc';
if mod(size(x,2),2)
    colors=colors_idx(1:size(x,2));
else
    colors=[colors_idx(1:size(x,2)/2) 'g' colors_idx(size(x,2)/2+1:size(x,2)) repmat('g',1,size(x,2)+1)];
end
if size(groups,1)>1; groups=groups'; end
separator=max(diff(sort(groups)))/size(x,2)/2;
h=figure;
hold on;
new_groups=[];
new_x=[];
new_labels=[];
for i=1:size(x,2)
    new_x=[new_x; x(:,i)];
    new_groups=[new_groups groups+(i-1)*separator];
    tmp=cell(size(x,1),1);
    tmp(:)={''};
    new_labels=[new_labels;tmp];
end
idx=unique(groups);
for i=1:length(idx)
    new_x=[new_x; NaN(size(x,2),1)];
    new_groups=[new_groups arrayfun(@(x) idx(i)+x*separator, size(x,2):size(x,2)*2-1)];
    tmp=cell(size(x,2),1);
    tmp(:)={''};
    new_labels=[new_labels;tmp];
    if ~mod(size(x,2),2)
        new_x=[new_x; NaN; NaN];
        new_groups=[new_groups idx(i)+separator*(size(x,2)-1)/2 idx(i)+separator*(size(x,2)/2*5+size(x,2)/2-1)/2];
        new_labels=[new_labels;g_labels{i};{''}];
    else
        new_labels(new_groups==(idx(i)+separator*floor(size(x,2)/2)))=g_labels(i);
    end
end
boxplot(new_x,new_groups,'colors',colors,'symbol','k.','labels',new_labels);
set(gca,'ticklength',[0 0]);

if nargin>4 % only for pairwise comparisons now
    sig_groups=zeros(length(unique(groups)),size(x,2));
    for i=1:size(sig_groups,1)
        sig_groups(i,:)=[1:size(x,2)/2 size(x,2)/2+2:size(x,2)+1]+(size(x,2)*2+2)*(i-1);
    end
    sig_groups=mat2cell(sig_groups,ones(1,length(unique(groups))),size(x,2));
    sigstar(sig_groups,pval);
end

if nargin>3
    for i=1:size(x,2)
        dummies(i)=plot([NaN NaN],[colors_idx(i) 's'],'markerfacecolor',colors_idx(i));
    end
    legend(dummies,v_labels,'fontsize',12,'location','southeast');
%     legend('boxoff');
end