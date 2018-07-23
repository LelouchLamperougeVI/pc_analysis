function [C,pval]=template_match(template,match,shuffles)
n=size(match,2);

idx=repmat((1:size(template,1))',1,size(match,1)-size(template,1)+1);
idx=idx+(0:(size(match,1)-size(template,1)));
match=match';
temp=match(:,idx);
match=reshape(temp,size(template,1)*size(match,1),size(match,2)-size(template,1)+1);

template=template';
template=template(:);

[C,pval]=fast_spearman(match,template,n,shuffles,false);