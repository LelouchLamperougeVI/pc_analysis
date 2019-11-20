function enregistrer(obj, overwrite)
if nargin < 2
    overwrite = false;
end

lfp = obj;
fn = fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']);
if ~exist(fn, 'file') || overwrite
    save(fn, 'lfp');
end

analysis = obj.analysis;
fn = fullfile(obj.session.wd, ['analysis' num2str(obj.session.id) '_plane' num2str(obj.twop.plane) '.mat']);
if ~isempty(obj.analysis) && (~exist(fn, 'file') || overwrite)
    save(fn, 'analysis');
end