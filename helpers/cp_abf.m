function animals = cp_abf(varargin)
% Copy remote abf files to local analysis folder
% {'scr', 'source', 'remote'}       dir     - the remote location
% {'dst', 'destination', 'local'}   dir     - the local analysis folder

ops = parse_inputs(varargin);

isdir = 1;
regex = {'^\d\d\d\d_\d\d_\d\d$'};
fcn = {@(x) datetime(x, 'inputformat','yyyy_MM_dd')};
animals = dirTree(ops.local, 'isdir',isdir, 'regex',regex, 'fcn',fcn);
while ~xor(isempty(animals.children), isempty(animals.var))
    isdir = [NaN, isdir];
    regex = [NaN, regex];
    fcn = [NaN, fcn];
    animals = dirTree(ops.local, 'isdir',isdir, 'regex',regex, 'fcn',fcn);
end

paths = animals.fullfile;

% for ii = 1:length(paths)
%     abf = fullfile(ops.remote, animals(ii).name, [animals(ii).date '_*.abf']);
%     try
%         copyfile(abf, animals(ii).parent);
%         disp(['copied: ' abf])
%     catch
%         warning(['failed to copy: ' abf]);
%     end
% end

function ops = parse_inputs(inputs)

ops.local = pwd;
ops.remote = '/mnt/mohaj1/homes/dun.mao/behavior';

count = 1;
while count <= length(inputs)
    switch inputs{count}
        case {'src', 'source', 'remote'}
            ops.remote = inputs{count + 1};
        case {'dst', 'destination', 'local'}
            ops.local = inputs{count + 1};
        otherwise
            error(['''' inputs{count} ''' is not a valid parameter']);
    end
    count = count + 2;
end