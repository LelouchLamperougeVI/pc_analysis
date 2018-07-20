function C=template_match(template,match)

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

% C=corr(match,template,'type','Spearman');
C=fast_spearman(match,template,false);