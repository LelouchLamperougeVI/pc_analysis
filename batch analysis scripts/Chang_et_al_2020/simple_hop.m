function retrieved = simple_hop(patterns, partials)
% Implementation of a simple synchronous, symmetric, discrete and
% fully-connected Hopfield network with Hebbian learning rule.

% patterns = ...
%    [1,-1,-1,-1,1,1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,1,-1,-1,-1,1
%    -1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,-1,1,1
%    1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,-1,-1,-1
%    -1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1];
% patterns = patterns';
% patterns = -patterns;
% 
% partials = ...
%     [1,-1,-1,-1,-1,1,1,1,1,1,-1,1,1,1,1,-1,1,1,1,1,1,1,-1,-1,-1];
% partials = partials';
% partials = -partials;

max_iter = 50;

w = patterns .* permute(patterns, [3 2 1]);
w = sum(w, 2);
w = w ./ size(patterns, 2);
w = squeeze(w);
w = w .* ~diag(ones(size(w, 1), 1));

retrieved = nan(size(partials, 2), 1);
for ii = 1:size(partials, 2)
    s = partials(:, ii);
    
    count = 1;
    while all(sum(patterns ~= s) ~= 0)
        s = iter(w, s);
        [~, retrieved(ii)] = min(sum(patterns ~= s));
        count = count + 1;
        if count > max_iter
            retrieved(ii) = nan;
            break
        end
    end
end

function s = iter(w, s)
s = double(s);
s = w * s;
s = double(s > 0);
s(~s) = -1;