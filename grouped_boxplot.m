function h = grouped_boxplot(x,groups,g_labels,v_labels)
% Grouped box plot with multiple variables
% Inputs
%   x: matrix with R observations and C variables
%   groups: grouping numeric vector
%   g_labels: group labels
%   v_labels: variable labels

colors_idx='rbkmc';
colors=colors_idx(1:size(x,2)+1);
if size(groups,1)>1; groups=groups'; end
separator=max(diff(sort(groups)))/size(x,2)/2;
h=figure;
new_groups=[];
new_x=[];
new_labels=[];
for i=1:size(x,2)
    new_x=[new_x; x(:,i)];
    new_groups=[new_groups groups+(i-1)*separator];
end
idx=unique(groups);
for i=1:length(idx)
    new_x=[new_x; NaN(size(x,2),1)];
    new_groups=[new_groups arrayfun(@(x) idx(i)+x*separator, size(x,2):size(x,2)*2-1)];
    if ~mod(size(x,2),2)
        new_x=[new_x; NaN; NaN];
        new_groups=[new_groups idx(i)+separator*size(x,2)/4 idx(i)+separator*size(x,2)*3/4];
    end
    tmp=cell(floor(size(x,2)/2),1);
    tmp{:}={''};
    new_labels=[new_labels tmp g_labels{i} tmp];
end
boxplot(new_x,new_groups,'colors',colors,'symbol','k.');