function I = dendro_colours(Z, clust)
% find the indices of the cluster in the dendrogram

I = [];
if(length(clust)>1)
    pairs = ones(length(clust));
    for i = clust
        for j = clust
            d = all(Z(:,1:2) == [i j], 2);
            if sum(d)
                pairs(clust==i,clust==j) = Z(d,3);
            end
        end
    end
    
    [d,idx] = min(pairs(:));
    d = find(Z(:,3)==d);
    [i,j] = ind2sub(size(pairs), idx);
    clust([i j]) = [];
    clust = [clust d+size(Z,1)+1];
    
    I = [I d];
    I = [I dendro_colours(Z, clust)];
end