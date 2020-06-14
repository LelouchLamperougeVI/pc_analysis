function rm_redund(obj)
% Remove overlapping/redundant neurons in multi-plane recordings

neur_diam = 30; % assuming that the diameter of a neuron is no more than this (this estimate should aim for the lower bound)

if ~isfield(obj.topo, 'centroid')
    warning('Topography hasn''t been build yet. Attempting to run obj.topography()');
    obj.topography;
end

% only across planes
d = abs(obj.twop.planes.depth - obj.twop.planes.depth');
d = d < obj.twop.planes.maxStep & d > 0;

% correlation between time-series
deconv = obj.twop.deconv;
deconv = deconv(~isnan(deconv));
deconv = reshape(deconv, [], size(obj.twop.deconv,2));
r = corr(deconv);

d = d & obj.topo.distances < neur_diam & r > obj.twop.planes.maxR;

d = triu(d, 1);

while sum(d, 'all') > 0
    [r, t] = find(d, 1);
    O = sum(any(obj.topo.maskNeurons == r, 3) .* any(obj.topo.maskNeurons == t, 3) , 'all') /...
        ( sum(any(obj.topo.maskNeurons == r, 3) , 'all') + sum(any(obj.topo.maskNeurons == t, 3) , 'all') - sum(any(obj.topo.maskNeurons == r, 3) .* any(obj.topo.maskNeurons == t, 3) , 'all') );
    if O > obj.twop.planes.ol
        if ( mean(deconv(:, r)) / std(deconv(:, r)) ) > ( mean(deconv(:, t)) / std(deconv(:, t)) )
            obj.rm_neurons(t);
            d(t, :) = [];
            d(:, t) = [];
        else
            obj.rm_neurons(r);
            d(r, :) = [];
            d(:, r) = [];
        end
    else
        d(r, t) = 0;
    end
end

