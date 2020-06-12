function rm_neurons(obj, n)
% Safely remove neurons

obj.twop.deconv(:, n) = [];
obj.twop.planes.plane_members(n) = [];
obj.topo.maskNeurons( ismember(obj.topo.maskNeurons, n) ) = 0;

try
    obj.topo.centroid(:, n) = [];
    obj.topo.distances(:, n) = [];
    obj.topo.distances(n, :) = [];
catch
end