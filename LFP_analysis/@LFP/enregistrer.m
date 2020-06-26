function enregistrer(obj, overwrite, transferable)
% Save the analysis structure
% To overwrite existing analysis, obj.enregistrer(true)
% To save as generic analysis file shareable across session,
% obj.enregistrer(false, true)
% ^typically used for RRR

% DEPRICATED

if nargin < 2
    overwrite = false;
end
if nargin < 3
    transferable = false;
end

% lfp = obj;
% fn = fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']);
% if ~exist(fn, 'file') || overwrite
%     save(fn, 'lfp');
% end

analysis = obj.analysis;
if transferable
    fn = fullfile(obj.session.wd, 'analysis.mat');
else
    fn = fullfile(obj.session.wd, ['analysis' num2str(obj.session.id) '_plane' strjoin(strsplit(num2str(obj.twop.planes.planes)), '_') '.mat']);
end
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