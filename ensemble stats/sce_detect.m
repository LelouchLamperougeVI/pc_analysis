function sce=sce_detect(assemblies,deconv)
shuffles=100;

sce=zeros(size(deconv,1),length(assemblies));
for i=1:length(assemblies)
    A=deconv(:,assemblies{i});
    [p,~,occ]=unique(A,'rows');
    occ=accumarray(occ,1);
    
    s_occ=zeros(length(occ),shuffles);
    for j=1:shuffles
        shuffled=mat_circshift(A,randi(size(A,1),1,size(A,2)));
        s_occ(:,j)=arrayfun(@(x) sum(all(p(x,:)==shuffled,2)),1:size(p,1));
    end
    
    occ=occ./sum(occ);
    s_occ=s_occ./sum(s_occ);
    
    idx=occ<prctile(s_occ',99)';
    idx=p(idx,:);
    idx=arrayfun(@(x) all(idx(x,:)==A,2), 1:size(idx,1), 'uniformoutput',false);
    idx=cell2mat(idx);
    sce(:,i)=~any(idx,2);
end
