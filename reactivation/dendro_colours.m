function I = dendro_colours(Z, clust)
% find the indices of the cluster in the dendrogram
% TODO: you and I both know that I hate recursions are they are slow as fuck...
% get rid of it! or at least the for-loops

I = [];
if(length(clust)>1)
%     pairs = ones(length(clust));
%     for i = clust
%         for j = clust
%             d = all(Z(:,1:2) == [i j], 2);
%             if sum(d)
%                 pairs(clust==i,clust==j) = Z(d,3);
%             end
%         end
%     end
    [~, i] = intersect(Z(:,1), clust);
    [~, j] = intersect(Z(:,2), clust);
    idx = intersect(i,j);
    
    mask = nan(size(Z,1),1);
    mask(idx) = 1;
    
    [~,idx] = min(Z(:,3) .* mask);
    
%     [d,idx] = min(pairs(:));
%     d = find(Z(:,3)==d);
%     [i,j] = ind2sub(size(pairs), idx);
%     clust([i j]) = [];
%     clust = [clust d+size(Z,1)+1];

    clust(clust==Z(idx,1)) = [];
    clust(clust==Z(idx,2)) = [];
    clust = [clust idx+size(Z,1)+1];
    
%     I = [I d];
    I = [I idx];
    I = [I dendro_colours(Z, clust)];
end