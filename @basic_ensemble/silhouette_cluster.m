function [clust, s, ticks] = silhouette_cluster(obj, Z, D, e_size)

if nargin < 2
    Z = obj.tree;
    D=1-abs(obj.R);
    D(1:length(D)+1:numel(D))=0;
    
    e_size = obj.ops.e_size;
end

ticks=[];
last_clust = 1;
last_max = 0;
s=[];
for i = 2:length(D)
    c = cluster(Z, 'maxclust', i);
    A = accumarray(c, 1);
    A = A > e_size;
    if sum(A) >= last_clust
        last_clust = sum(A);
        A = A .* (1:length(A))';
        A(~A) = [];
        A = ismember(c, A);
        s = [s mean(silhouette([], c(A), squareform(D(A,A))))];
        ticks = [ticks last_clust];
        if s(end) > last_max
            last_max = s(end);
            clust = c;
        end
    end
end