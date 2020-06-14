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

for ii = 1:length(ass.clust)
    plot3(ass.topo.centroid(1, ass.clust{ii}), ass.topo.centroid(2, ass.clust{ii}), -ass.twop.planes.depth(ass.clust{ii}), '.', 'markerfacecolor', ass.colours(ii, :), 'markersize', 30);
end

xlabel('galvo')
ylabel('resonant')
zlabel('piezo')


