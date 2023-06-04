function md = ccpc_ll(obj)

behavior = obj.analysis.behavior;
ops.unit_pos=behavior.unit_pos;
ops.unit_vel=behavior.unit_vel;
ops.frame_ts=behavior.frame_ts;
ops.trials=behavior.trials;

thres=noRun(ops.unit_vel);

thres=(ops.unit_vel>thres | ops.unit_vel<-thres) & (ops.trials(1) < ops.frame_ts & ops.trials(end) > ops.frame_ts);
ops.unit_pos=ops.unit_pos(thres);
ops.frame_ts=ops.frame_ts(thres);

deconv = obj.analysis.original_deconv;
deconv = ca_filt(deconv);
n = deconv(thres,:);
x = ops.unit_pos(:);
x_ = arrayfun(@(ii) x(n(:, ii) > 0), 1:size(n, 2), 'UniformOutput', false);
n_ = arrayfun(@(ii) n(n(:, ii) > 0, ii), 1:size(n, 2), 'UniformOutput', false);

belt;
blt = linspace(min(x), max(x), length(belt_num) + 1);
blt = discretize(x, blt);
blt = belt_num(blt);
c_ = arrayfun(@(ii) blt(n(:, ii) > 0), 1:size(n, 2), 'UniformOutput', false);
compartments = unique(belt_num);
c_ = cellfun(@(x) categorical(x, compartments), c_, 'UniformOutput', false);
compartments = categorical(compartments);
c_ = cellfun(@double, c_, 'UniformOutput', false);

md.ll = zeros(length(x_), 3);
md.lambda = zeros(length(x_), length(compartments));
md.bad = false(length(x_), 1);

% % without baseline
% md.pc = zeros(length(x_), 3);
% gaussfun = @(v, x) v(1) .* exp(-(x - v(2)).^2 ./ (2 * v(3)^2));
% A = [-1  0  0
%       0  1  0
%       0 -1  0
%       0  0 -1];
% bnd = [0
%        max(x)
%        min(x)
%        0];
% with baseline
md.pc = zeros(length(x_), 4);
gaussfun = @(v, x) v(1) .* exp(-(x - v(2)).^2 ./ (2 * v(3)^2)) + v(4);
A = [-1  0  0  0
      0  1  0  0
      0 -1  0  0
      0  0  1  0
      0  0 -1  0
      0  0  0 -1
     -1  0  0  1]; % peak amplitude higher than baseline
bnd = [0
       max(x)
       min(x)
       range(x)/4
       0
       0
       0];

for ii = 1:length(x_)
    % gaussian model
    x = x_{ii}(:);
    n = n_{ii}(:);
%     x0 = [mean(n), sum(n.*x)/sum(n), 10]; % without baseline
    x0 = [mean(n), sum(n.*x)/sum(n), 10, median(n)]; % with baseline
    fun = @(v) - sum( n .* log(gaussfun(v, x)) - gaussfun(v, x) - gammaln(n + 1) );
    md.pc(ii, :) = fmincon(fun, x0, A, bnd, [], [], [], [], [], optimoptions('fmincon', 'display', 'off'));
    lambda = gaussfun(md.pc(ii, :), x);
    md.ll(ii, 1) = sum(n .* log(lambda) - lambda - gammaln(n + 1));
    
    % cue compartment model
    lambda = accumarray(c_{ii}(:), n_{ii}(:), [length(compartments), 1], @mean);
    md.ll(ii, 2) = sum(n_{ii}(:) .* log(lambda(c_{ii})) - lambda(c_{ii}) - gammaln(n_{ii}(:) + 1));
    md.lambda(ii, :) = lambda;
    
    % null model
    lambda = mean(n_{ii});
    md.ll(ii, 3) = sum(n_{ii}(:) .* log(lambda) - lambda - gammaln(n_{ii}(:) + 1));
end

md.ratio = -2 .* (md.ll(:, 3) - md.ll(:, 1:2));
md.ssize = cellfun(@length, n_);