function len = ctnseg(a)
% Find length of longest consecutive segment in logical array

if ~ismatrix(a)
    error('So sorry, can''t do more than 2 dimensions');
end
if size(a,1) == 1
    a = a';
end

len = zeros(size(a,2), 1);
seg = diff(a);
for ii = 1:size(a,2)
    idx = find(seg(:, ii));
    if isempty(idx)
        continue
    end
    if seg(idx(1), ii) ~= 1
        idx = [1; idx];
    end
    if seg(idx(end), ii) ~= -1
        idx = [idx; size(seg, 1)];
    end
    len(ii) = max(diff(idx) + 1);
end