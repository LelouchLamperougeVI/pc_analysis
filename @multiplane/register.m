function overlap = register(obj, a, b, varargin)
% Find overlapping and differing ROIs between 2 sessions
% 'Name, Value' pairs:
%   'overlap' / 'threshold'     fraction of overlap - .5 (default)
%   'score' / 'method'          method for calculating overlap score - 'one-sided' / 'one' or 'two-sided' / 'two' (default)
%   'plot' / 'plotFlag',        false (default) - true
%
% one-sided: $O = frac{ overlap pixels }{ moving ROI pixels }$
% two-sided: $O = frac{ overlap pixels }{ ROI1 pixels + ROI2 pixels - overlap pixels}$
%
% Output: the ROIs of 'b' that matchs to 'a'

% check if same day/session
if ( obj.planes(a).session.date == obj.planes(b).session.date ) && ( size(obj.planes(a).twop.deconv, 2) == size(obj.planes(b).twop.deconv, 2) )
    disp(['Planes ' num2str(a) ' and ' num2str(b) ' were acquired on the same day']);
    overlap.roi = num2cell( 1:size(obj.planes(a).twop.deconv, 2) );
    overlap.score = num2cell( ones(1, size(obj.planes(a).twop.deconv, 2)) );
    return;
end

ops = parse_inputs(varargin);

maskA = obj.planes(a).topo.maskNeurons;
maskB = obj.planes(b).topo.maskNeurons;

fixed = double(~~maskA);
moving = double(~~maskB);

[optimizer, metric] = imregconfig('monomodal');
transform = imregtform(moving, fixed, 'rigid', optimizer, metric);
maskReg = imwarp(maskB, transform, 'nearest', 'outputview', imref2d(size(fixed)), 'fillvalues', 0);

stack = maskReg == permute(1:max(maskReg(:)), [1 3 2]); % tried to be lazy, ended up overcomplicating things...
stack = maskA .* stack;
stack = reshape(stack, [numel(maskA) size(stack,3)]);
stack = arrayfun(@(x) accumarray(stack(:, x)+1, 1, [max(maskA(:))+1 1]), 1:size(stack,2), 'uniformoutput',false);
stack = cell2mat(stack);
stack = stack(2:end,:);

movPix = accumarray(maskReg(:)+1, 1, [max(maskReg(:))+1 1]);
movPix = movPix(2:end);
if strcmp(ops.method, 'one')
    stack = stack ./ movPix';
    overlap.roi = arrayfun(@(x) find(stack(:,x) > ops.threshold), 1:size(stack,2), 'uniformoutput',false);
    overlap.score = arrayfun(@(x) stack(overlap.roi{x}, x), 1:size(stack,2), 'uniformoutput',false);
elseif strcmp(ops.method, 'two')
    fixPix = accumarray(maskA(:)+1, 1, [max(maskA(:))+1 1]); % number of pixels for fixed ROIs
    fixPix = fixPix(2:end);
    
    stack = stack ./ ( fixPix + movPix' - stack );
    overlap.roi = arrayfun(@(x) find(stack(:,x) > ops.threshold), 1:size(stack,2), 'uniformoutput',false);
    overlap.score = arrayfun(@(x) stack(overlap.roi{x}, x), 1:size(stack,2), 'uniformoutput',false);
end

overlap.ops = ops;

disp(['Detected ' num2str(sum(~cellfun(@isempty, overlap.roi))) '/' num2str(length(overlap.roi)) ' stable neurons']);

if ops.plot
    registered = imwarp(moving, transform, 'nearest', 'outputview', imref2d(size(fixed)), 'fillvalues', 0);
    figure
    ax(1) = subplot(1,3,1);
    imshowpair(fixed, moving,'Scaling','joint')
    title('original')
    axis xy
    axis square
    ax(2) = subplot(1,3,2);
    imshowpair(fixed, registered,'Scaling','joint')
    title('registered')
    axis xy
    axis square
    
    idx = find(~cellfun(@isempty, overlap.roi));
    idx = ismember(maskB, idx);
    [y,x] = find(idx);
    ax(3) = subplot(1,3,3);
    plot(x,y,'r.');
    hold on
    idx = find(cellfun(@isempty, overlap.roi));
    idx = ismember(maskB, idx);
    [y,x] = find(idx);
    plot(x,y,'k.');
    legend('overlap', 'differ')
    axis square
    
    linkaxes(ax)
end


function ops = parse_inputs(inputs)
ops.plot = false;
ops.method = 'two';
ops.threshold = .5;

count = 1;
while(count < length(inputs))
    switch lower(inputs{count})
        case {'plot', 'plotflag'}
            ops.plot = inputs{count+1};
        case {'score', 'method'}
            switch lower(inputs{count+1})
                case {'one', 'one-sided'}
                    ops.method = 'one';
                case {'two', 'two-sided'}
                    ops.method = 'two';
            end
        case {'overlap', 'threshold'}
            ops.threshold = inputs{count+1};
        otherwise
            error([inputs{count} ' is not a valid ''name'' specifier']);
    end
    
    count = count + 2;
end