function enregistrer(obj,fn)
lfp = obj;
if nargin<2
    fn = fullfile(obj.session.wd, ['lfp' num2str(obj.session.id) '.mat']);
end
save(fn, 'lfp');