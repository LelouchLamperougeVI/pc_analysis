function save(obj, target, spec)

if nargin < 3
    spec = false;
end

switch lower(target)
    case 'analysis'
        analysis = obj.analysis;
        if spec
            fn = fullfile(obj.session.wd, ['analysis' num2str(obj.session.id) '_plane' strjoin(strsplit(num2str(obj.twop.planes.planes)), '_') '.mat']);
        else
            fn = fullfile(obj.session.wd, 'analysis.mat');
        end
        if isempty(analysis)
            error(['Analysis is missing. Please run ' class(obj) '.perform_analysis()']);
        end
        if exist(fn, 'file')
            warning('Existing analysis file will be overwritten.');
        end
        save(fn, 'analysis');
        
    case {'chan', 'channel', 'channels'}
        chan = obj.abf.Channels;
        if spec
            fn = fullfile(obj.session.wd, ['channels' num2str(obj.session.id) '_plane' strjoin(strsplit(num2str(obj.twop.planes.planes)), '_') '.mat']);
        else
            fn = fullfile(obj.session.wd, 'channels.mat');
        end
        if exist(fn, 'file')
            warning('Existing analysis file will be overwritten.');
        end
        save(fn, 'chan');
        
    otherwise
        error([target 'is undefined.']);
end