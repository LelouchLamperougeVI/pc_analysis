function plot(obj, type)

switch(lower(type))
    case 'masks'
        figure;
        k = 1;
        for ii = 1 : size(obj.ROI.overlap, 1)
            for jj = 1:size(obj.ROI.overlap, 2)
                subplot(size(obj.ROI.overlap, 1), size(obj.ROI.overlap, 2), k);
                imshowpair(~~obj.planes(ii).topo.maskNeurons, ~~obj.ROI.overlap{ii, jj}.maskReg);
                title({['date: ' datestr(obj.planes(ii).session.date)], ['X' datestr(obj.planes(jj).session.date)]});
                k = k + 1;
            end
        end
end