clusts = 1:length(ass.clust);
clusts = 15;

figure
hold on

x = linspace(0, ass.topo.FOV(1), size(ass.topo.mimg, 1));
y = linspace(0, ass.topo.FOV(2), size(ass.topo.mimg, 2));
sections = -unique(ass.twop.planes.depth);
[x, y, z] = meshgrid(y, x, sections);
slice(x, y, z, ass.topo.mimg, [], [], sections)
axis image
colormap gray
shading interp
alpha(.3)

for ii = clusts
    plot3(ass.topo.centroid(1, ass.clust{ii}), ass.topo.centroid(2, ass.clust{ii}), -ass.twop.planes.depth(ass.clust{ii}), '.', 'markerfacecolor', ass.colours(ii, :), 'markersize', 30);
end

xlabel('galvo')
ylabel('resonant')
zlabel('piezo')


%%
figure
subplot(2, 2, 1);
cdfplot(cell2mat(obj.topo.clust.zsilhouette));
xlim([-1 1]);
xlabel('silhouette score')
ylabel('cum. freq.')
xline(0);

subplot(2, 2, 2);
violin(obj.topo.clust.zsilhouette');
ylim([-1 1]);
xlabel('ensemble');
ylabel('silhouette score')
yline(0);

subplot(2, 2, 3);
d = cellfun(@(x) squareform(obj.topo.zdistances(x, x))', obj.clust, 'uniformoutput', false);
violin(d);
xlabel('ensemble');
ylabel('distance (\mum)')
d = squareform(obj.topo.zdistances);
yline(median(d));
title('line: sample median')

subplot(2, 2, 4);
depth = obj.twop.planes.depth(cell2mat(obj.clust));
s = cell2mat(obj.topo.clust.zsilhouette);
s = arrayfun(@(x) s(depth == x), unique(depth), 'uniformoutput', false);
violin(s, 'xlabel', strsplit(num2str(unique(depth))));
ylim([-1 1]);
xlabel('plane depth');
ylabel('silhouette score')
yline(0);


%%
clusts = 1:length(obj.clust);
% clusts = 15;
k=2;
for plane = 1:size(obj.topo.mimg, 3)
    if ~mod(plane-1, k^2)
        figure;
    end
    subplot(k, k, mod(plane-1, k^2) + 1);
    
    mimg = obj.topo.mimg(:, :, plane);
    mimg = (mimg - min(mimg(:))) ./ range(mimg(:));
    imshow(mimg)
    hold on

    mask = zeros(size(obj.topo.maskNeurons, 1), size(obj.topo.maskNeurons, 2), 3);
    for ii = clusts
        idx = ismember(obj.topo.maskNeurons(:, :, plane), obj.clust{ii});
        idx = repmat(idx, 1, 1, 3);
        mask = mask + double(idx) .* permute(obj.colours(ii, :), [1 3 2]);
    end

    h = imshow(mask);
    set(h, 'alphadata', .5);
    title(obj.twop.planes.plane_names{plane});
end





