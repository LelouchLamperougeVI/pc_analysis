function obj = set(obj, ops)

if isempty(obj.ops)
    obj.ops.discard = [0 1 12 13 14]; % planes to discard
end

count = 1;
while count < length(ops)
    switch lower(ops{count})
        case {'discard', 'garbage'}
            obj.ops.discard = ops{count + 1};
        otherwise
            error(['The argument ' ops{count} ' is not a valid property.']);
    end
    count = count + 2;
end