function enregistrer(obj, overwrite)
% save the analysis structure

if nargin < 2
    overwrite = false;
end

% lfp = obj;
% fn = fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']);
% if ~exist(fn, 'file') || overwrite
%     save(fn, 'lfp');
% end

analysis = obj.analysis;
fn = fullfile(obj.session.wd, ['analysis' num2str(obj.session.id) '_plane' strjoin(strsplit(num2str(obj.twop.planes.planes)), '_') '.mat']);
% if ~isempty(obj.analysis) && (~exist(fn, 'file') || overwrite)
%     save(fn, 'analysis');
% end
if isempty(analysis)
    error(['Analysis is missing. Please run ' class(obj) '.perform_analysis()']);
end

if exist(fn, 'file') && overwrite
    warning('Existing analysis file will be overwritten.');
elseif exist(fn, 'file') && ~overwrite
    error('Analysis file already exists. Please confirm overwrite.');
end

save(fn, 'analysis');