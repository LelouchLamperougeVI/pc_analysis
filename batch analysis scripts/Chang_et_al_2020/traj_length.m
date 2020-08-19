function [l, s, e] = traj_length(trajectories, bl_ratio)
% compute the length, start and end of trajectories

if nargin < 2
    bl_ratio = 150 / 50;
end

l = cell(size(trajectories, 2), 1);
s = cell(size(trajectories, 2), 1);
e = cell(size(trajectories, 2), 1);
% idx = diff(trajectories);
idx = diff(repmat(trajectories, 3, 1));
idx = idx(size(trajectories, 1):size(trajectories, 1)*2-1, :);

for ii = 1:size(trajectories, 2)
    starts = find(idx(:, ii) == 1);
    ends = find(idx(:, ii) == -1);
    if idx(find(idx(:, ii) ~= 0, 1), ii) == -1
        starts = circshift(starts, 1);
    end
    lengths = min(mod([starts-ends ends-starts], size(trajectories, 1)), [], 2);
    
    l{ii} = lengths;
    s{ii} = starts;
    e{ii} = ends;
end

% for ii = 1:size(trajectories, 2)
%     if ~sum(trajectories(:, ii))
%         continue
%     end
%     if idx(find(idx(:, ii) ~= 0, 1), ii) == -1
%         l{ii} = cat(1, l{ii}, find(idx(:, ii) == -1, 1));
%         s{ii} = cat(1, s{ii}, 0);
%         e{ii} = cat(1, e{ii}, find(idx(:, ii) == -1, 1));
%         idx(find(idx(:, ii) == -1, 1), ii) = 0;
%     end
%     if idx(find(idx(:, ii) ~= 0, 1, 'last'), ii) == 1
%         l{ii} = cat(1, l{ii}, size(trajectories, 1) - find(idx(:, ii) == 1, 1, 'last'));
%         s{ii} = cat(1, s{ii}, find(idx(:, ii) == 1, 1, 'last'));
%         e{ii} = cat(1, e{ii}, size(trajectories, 1));
%         idx(find(idx(:, ii) == 1, 1, 'last'), ii) = 0;
%     end
%     l{ii} = cat(1, l{ii}, find(idx(:, ii) == -1) - find(idx(:, ii) == 1));
%     s{ii} = cat(1, s{ii}, find(idx(:, ii) == 1));
%     e{ii} = cat(1, e{ii}, find(idx(:, ii) == -1));
%     
%     l{ii} = l{ii} .* bl_ratio;
% end