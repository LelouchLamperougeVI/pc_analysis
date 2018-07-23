function [C,pval]=template_match(template,match,shuffles)

s=1:size(template,1);
idx=[];
for i=1:size(match,1)-size(template,1)
    idx=[idx s+i];
end
match=match';
temp=match(:,idx);
match=reshape(temp,size(template,1)*size(match,1),size(match,2)-size(template,1));

template=template';
template=template(:);

[C,pval]=fast_spearman(match,template,shuffles,false);