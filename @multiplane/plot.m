function plot(obj, type)

switch(lower(type))
    case 'masks'
        figure;
        k = 1;
        for ii = obj.ops.target.index
            for jj = obj.ops.target.index
                subplot(length(obj.ops.target.index), length(obj.ops.target.index), k);
                imshowpair(~~obj.planes(ii).topo.maskNeurons, ~~obj.ROI.overlap{ii, jj}.maskReg);
                title({['date: ' datestr(obj.planes(ii).session.date)], ['X' datestr(obj.planes(jj).session.date)]});
                k = k + 1;
            end
        end
end