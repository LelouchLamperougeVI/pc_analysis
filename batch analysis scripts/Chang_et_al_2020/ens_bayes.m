function md = ens_bayes(obj)
% bayesian decoding of cue and position by ensemble
belt;

behavior = obj.analysis.behavior;
ops.unit_pos=behavior.unit_pos;
ops.unit_vel=behavior.unit_vel;
ops.frame_ts=behavior.frame_ts;
ops.trials=behavior.trials;

thres=noRun(ops.unit_vel);
thres=(ops.unit_vel>thres | ops.unit_vel<-thres) & (ops.trials(1) < ops.frame_ts & ops.trials(end) > ops.frame_ts);

idx = knnsearch(ops.frame_ts(:), ops.trials(:));
idx = cat(1, 1, idx(:)); idx(end + 1) = length(ops.unit_pos) + 1;
idx = repelem((0:length(idx) - 2)', diff(idx));

ops.unit_pos=ops.unit_pos(thres);
trials = idx(thres);
ops.frame_ts=ops.frame_ts(thres);

deconv = obj.analysis.original_deconv;
deconv = ca_filt(deconv);
n = deconv(thres,:);
x = ops.unit_pos(:);
trials = trials(:);

blt = linspace(min(x), max(x), length(belt_num) + 1);
blt = discretize(x, blt);
blt = belt_num(blt);

bins = 150;
clusts = obj.ensembles.clust;
err = zeros(bins, length(clusts));
acc = zeros(1, length(clusts));
individual = zeros(length(unique(blt)) - 1, length(clusts));
% true = zeros(length(x), length(clusts)); decoded = true;
for ii = 1:length(clusts)
    md = sam_mape(x, n(:, clusts{ii}), trials, 'dt', round(obj.twop.fs), 'bins', bins, 'circ', true, 'cv', true', 'mle', false, 'penalty', eps, 'smooth', 3);
%     true(:, ii) = md.x;
%     decoded(:, ii) = md.decoded;
    err(:, ii) = md.err;
    
    md = sam_mape(blt, n(:, clusts{ii}), trials, 'dt', round(obj.twop.fs), 'bins', length(unique(blt)), 'circ', false, 'cv', true', 'mle', false, 'penalty', eps, 'smooth', 0);
    temp = accumarray([double(categorical(md.x(:), unique(md.x))), double(categorical(md.decoded(:), unique(md.x)))], 1, [1, 1] .* length(unique(blt)), @sum);
    [~, idx] = max(sum(temp, 2));
    idx = setxor(1:length(temp), idx);
    acc(ii) = sum(diag(temp(idx, idx))) ./ sum(temp(idx, :), 'all');
    individual(:, ii) = diag(temp(idx, idx)) ./ sum(temp(idx, :), 2);
end

clear md
md.acc = acc;
md.individual = individual;
md.err = err;