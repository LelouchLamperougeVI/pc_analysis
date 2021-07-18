function md = pc_cc_simanneal(n, x, dt, varargin)
% Fit "place cells" and "cue compartment" models to neuronal data by
% simulated annealing.
%
% Inputs: pc_cc_simanneal(n, x, blt, dt, 'name', 'value')
%   n                   deconvolved calcium fluorescence (Time x Neurons)
%   x                   linear position vector of length Time
%   dt                  delta-t parameter for Poisson log-likelihood
%   'plot', false       plot model fits
%   'prog', true        display progress bar
%   'reject', vect      logical vector of rejection epochs of length Time
%   'bins', 50          number of bins for discretized position
%   'blt', blt          belt vector (default to Ingrid's belt)
%   'comp', 'intra'     compartmentalization model
%       'full':     both inter- and intra- cue spaces represent a discrete
%                   firing rate compartment
%       'intra':    each intra-cue space has a discrete firing rate, while
%                   the inter-cue space has a static baseline lambda
%       'bin':      old-school centre-on/surround-off model

% parse them name,value pair
ops = parse_ops(varargin);

% discretize position
edges = linspace(min(x), max(x), ops.bins + 1);
x = discretize(x, edges);
x = x(:);

% window summing firing rates (for log-likelihood estimation)
kernel = ones(dt, 1);
n = arrayfun(@(x) conv(n(:, x), kernel, 'same'), 1:size(n, 2), 'UniformOutput', false);
n = cell2mat(n);

% n = fast_smooth(n, dt);

% remove rejection epochs
n(ops.reject, :) = [];
x(ops.reject) = [];

% compartmentalize belt if not already
blt = double(~~ops.blt(:));
switch ops.comp
    case {'full'}
        blt = comparmentalize_full(blt);
    case {'intra'}
        blt = comparmentalize_intra(blt);
    case {'bin', 'binary'}
    otherwise
        error('Undefined compartmentalization model.');
end

% get stack
stack = arrayfun(@(ii) accumarray(x, n(:, ii), [ops.bins 1], @mean), 1:size(n, 2), 'UniformOutput', false);
stack = cellfun(@(x) x ./ dt, stack, 'UniformOutput', false);

% run model fitting
fit_pc = zeros(size(n, 2), 4);
fit_cc = zeros(size(n, 2), length(unique(blt)));

dq=parallel.pool.DataQueue;
if ops.prog
    h = waitbar(0,'model fitting... (simulated annealing)');
    afterEach(dq,@updateBar);
    progress = 0;
else
    h = [];
end
    function updateBar(~)
        waitbar(progress/size(n, 2), h);
        progress = progress + 1;
    end

parfor ii = 1:size(n, 2)
% for ii = 1:size(n, 2)
    n_single = n(:, ii);
    
    init = [max(stack{ii}),...
        find(stack{ii} == max(stack{ii}), 1),...
        sum(stack{1} > (mean(stack{1}) + std(stack{1})*2)),...
        mean(n_single)];
    max_t = [max(n_single) ops.bins ops.bins/4 mean(n_single)];
    init_t = [max(stack{ii})/5 ops.bins/20 init(3)/5 mean(n_single)];
%     options = optimoptions('simulannealbnd','InitialTemperature', init_t,...
%         'TemperatureFcn', 'temperaturefast', 'AnnealingFcn', 'annealingboltz',...
%         'ReannealInterval', 10, 'FunctionTolerance', 1e1,...
%         'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf, @saplottemperature, @saplotstopping});
    options = optimoptions('simulannealbnd','InitialTemperature', init_t,...
        'TemperatureFcn', 'temperatureboltz', 'AnnealingFcn', 'annealingboltz',...
        'ReannealInterval', 10, 'FunctionTolerance', 1e1);
    l = @(params) pc_md_opt_fun(params, x, n_single, dt);
    fit_pc(ii, :) = simulannealbnd(l, init, [0 0 0 0], max_t, options);
    
    init = accumarray(blt(x), n_single, [length(unique(blt)) 1], @mean);
    max_t = ones(length(unique(blt)), 1) .* max(n_single);
    init_t = ones(length(unique(blt)), 1) .* max(init) ./ 10;
    options = optimoptions('simulannealbnd','InitialTemperature', init_t,...
        'TemperatureFcn', 'temperatureboltz', 'AnnealingFcn', 'annealingboltz',...
        'ReannealInterval', 10, 'FunctionTolerance', 1e1);
%     options = optimoptions('simulannealbnd','InitialTemperature', init_t,...
%         'TemperatureFcn', 'temperatureboltz', 'AnnealingFcn', 'annealingboltz',...
%         'ReannealInterval', 10, 'FunctionTolerance', 1e1, ...
%         'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf, @saplottemperature, @saplotstopping});
    l = @(params) cc_md_opt_fun(params, x, n_single, dt, blt);
    fit_cc(ii, :) = simulannealbnd(l, init, zeros(length(unique(blt)), 1), max_t, options);
    
    send(dq, ii);
end
close(h)

% reconstruct neuronal activity with models
reconst_pc = arrayfun(@(ii) fit_pc(ii, 1) .* exp(-((1:(range(x)+1)) - fit_pc(ii, 2)).^2 ./ (2 * fit_pc(ii, 3)^2)) + fit_pc(ii, 4), 1:size(n, 2), 'UniformOutput', false); %reconst pos
reconst_pc = cell2mat(reconst_pc');
reconst_pc = arrayfun(@(ii) reconst_pc(ii, x), 1:size(n, 2), 'UniformOutput', false); % pos to time
reconst_pc = cell2mat(reconst_pc')';

reconst_cc = arrayfun(@(ii) fit_cc(ii, blt), 1:size(n, 2), 'UniformOutput', false); % belt to pos
reconst_cc = cell2mat(reconst_cc');
reconst_cc = arrayfun(@(ii) reconst_cc(ii, x), 1:size(n, 2), 'UniformOutput', false); % pos to time
reconst_cc = cell2mat(reconst_cc')';

% reconstruct stacks
md.stack_n = cell2mat(stack);
md.stack_pc = cell2mat(arrayfun(@(ii) accumarray(x, reconst_pc(:, ii), [ops.bins 1], @mean), 1:size(n, 2), 'UniformOutput', false));
md.stack_cc = cell2mat(arrayfun(@(ii) accumarray(x, reconst_cc(:, ii), [ops.bins 1], @mean), 1:size(n, 2), 'UniformOutput', false));

% (partial) explained variance
r_pc = corr(n, reconst_pc);
r_pc = diag(r_pc);
r_cc = corr(n, reconst_cc);
r_cc = diag(r_cc);
r = corr(reconst_pc, reconst_cc);
r = diag(r);

ev_pc = r_pc .^2;
ev_cc = r_cc .^2;

pev_pc = ( (r_pc - r_cc .* r) ./ (sqrt(1 - r_cc.^2) .* sqrt(1 - r.^2)) ) .^2;
pev_cc = ( (r_cc - r_pc .* r) ./ (sqrt(1 - r_pc.^2) .* sqrt(1 - r.^2)) ) .^2;

md.pc.params = fit_pc;
md.pc.reconst = reconst_pc;
md.pc.ev = ev_pc;
md.pc.pev = pev_pc;

md.cc.params = fit_cc;
md.cc.reconst = reconst_cc;
md.cc.ev = ev_cc;
md.cc.pev = pev_cc;

md.n = n;
md.x = x;
md.dt = dt;
md.ops = ops;

end

function ops = parse_ops(inputs)
ops.plot = false;
ops.prog = true;
ops.reject = [];
ops.bins = 50;
ops.comp = 'intra';
% Ingrid's belt object
ops.blt = [0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1,0,0];

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        case {'plot', 'plotflag'}
            ops.plot = logical(inputs{count+1});
        case {'prog', 'progress'}
            ops.prog = logical(inputs{count+1});
        case {'reject'}
            ops.reject = logical(inputs{count+1});
        case {'bins'}
            ops.bins = inputs{count+1};
        case {'blt', 'belt'}
            ops.blt = inputs{count+1};
        case {'comp', 'compartment', 'compartmentalization', 'comp_md'}
            ops.comp = inputs{count+1};
        otherwise
            error('undefined input sequence');
    end
    count = count + 2;
end
end

function l = pc_md_opt_fun(params, x, n, dt)
% loss function (negative log-likelihood) for place cells model
% l = -sum( n .* log(params(1)) - n .* (x - params(2)).^2 ./ (2*params(3)^2) - params(1) .*dt .* exp(-(x - params(2)).^2 ./ (2*params(3)^2)));

lambda = params(1) .* exp(-(x - params(2)).^2 ./ (2 * params(3) ^ 2)) + params(4);
l = -sum( n .* log(lambda) - dt .* lambda );
end

function l = cc_md_opt_fun(params, x, n, dt, blt)
% loss function (negative log-likelihood) for cue compartment model
x = blt(x);
x = x(:);
n = n(:);
params = params(:);

l = -sum( n .* log(params(x)) - dt .* params(x) );
end

function blt = comparmentalize_intra(blt)
% break belt into cue compartments
edges = [get_head(blt), get_head(blt(end:-1:1))];
edges = [find(edges(:, 1)), find(edges(end:-1:1, 2))];

temp = zeros(size(blt));
for ii = 1:length(edges)
    temp(edges(ii, 1) : edges(ii, 2)) = ii;
end
blt = temp + 1;
end

function blt = comparmentalize_full(blt)
% break belt into cue compartments
edges = find(abs(diff(blt)));
edges = cat(1, 0, edges(:), length(blt));

temp = zeros(size(blt));
for ii = 1:length(edges)-1
    temp((edges(ii) + 1) : edges(ii+1)) = ii;
end
blt = temp;
end