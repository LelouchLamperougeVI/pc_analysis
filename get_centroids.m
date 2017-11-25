function centroid_mask=get_centroids(analysis)
centroid_mask=zeros(size(analysis.maskNeurons));
for i=1:length(analysis.psth)
    centroids=regionprops(analysis.maskNeurons==i,'centroid');
    centroids=cat(1,centroids.Centroid);
    centroids=round(centroids);
    centroid_mask(centroids(1),centroids(2))=i;
end
centroid_mask=centroid_mask';