function assemblies=demix_overlaps(D,assemblies,thres)
% Search assemblies to find inter-clusters and intra-clusters overlaps
% between ensemble members

count=1;
for i=1:length(assemblies)
    Dm=D(assemblies{i},assemblies{i});
    Dm=squareform(Dm);
    Z=linkage(Dm,'weighted');
    idx=diff(Z(:,3))>(mean(diff(Z(:,3)))+3*std(diff(Z(:,3))));
    idx=find(idx);
    
    if isempty(idx)
        new_assemblies{count}=assemblies{i};
        count=count+1;
    else
        idx=[idx;size(Z,1)-1];
        temp=assemblies{i};
        for j=1:length(idx)
            clusters=cluster(Z,'criterion','distance','cutoff',Z(idx(j)+1,3));
            [new_assemblies{count},ia]=intersect(temp, assemblies{i}(clusters==mode(clusters)));
            temp(ia)=[];
            count=count+1;
        end
    end
end

assemblies=new_assemblies;

% count=1;
% for i=1:length(assemblies)
%     Dm=D(assemblies{i},assemblies{i});
%     Dm=squareform(Dm);
%     Z=linkage(Dm,'complete');
%     clusters=cluster(Z,'criterion','distance','cutoff',thres);
%     for j=1:max(clusters)
%         new_assemblies{count}=assemblies{i}(clusters==j);
%         count=count+1;
%     end
% end
% 
% assemblies=new_assemblies;

new_members=cell(1,length(assemblies));
for i=1:length(new_members); new_members{i}=[]; end

for i=1:length(assemblies)
    for j=i+1:length(assemblies)
        overlaps=D(assemblies{i},assemblies{j});
        overlaps=mean(overlaps)<thres & mean(overlaps,2)<thres;
        new_members{i}=[new_members{i} assemblies{j}(sum(overlaps)./size(overlaps,2) > .5)];
        new_members{j}=[new_members{j} assemblies{i}(sum(overlaps,2)./size(overlaps,1) > .5)];
    end
end
for i=1:length(assemblies)
    assemblies{i}=unique([assemblies{i} new_members{i}]);
end
        
fprintf('Detected assemblies: \n');
arrayfun(@(x) fprintf(['\t Assembly ' num2str(x) ':\t' mat2str(assemblies{x}) '\n']),1:length(assemblies));