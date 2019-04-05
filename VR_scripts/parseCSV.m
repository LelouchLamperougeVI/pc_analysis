function C = parseCSV(fn)
% All purpose parser for CSV files generated from TabletVR
% Assumes first value is timestamp

if nargin < 1
    [file,path] = uigetfile('*.csv');
    fn = fullfile(path,file);
end

full = fullfile(path, [file(1:end-4) '.mat']);
if exist(full, 'file')
    C = load(full);
    C = C.C;
    return;
end

fid = fopen(fn, 'r');

C = struct('temp','');

while ~feof(fid)
    s = fgetl(fid);
    s = split(s, ',');
    s{2} = regexprep(s{2}, '-', '');
    if ~isfield(C, s{2})
        C.(s{2}) = [];
    end
    idx=cellfun(@str2double, s)';
    if ~any(isnan(idx(3:end)))
        C.(s{2}) = [C.(s{2}); idx];
    else
        C.(s{2}) = [C.(s{2}); s'];
    end
end

C = rmfield(C, 'temp');

idx = fieldnames(C);
for i = 1:length(idx)
    C.(idx{i}) (:,2) = [];
end

C.date = datetime(file(3:16), 'inputformat', 'yyyyMMddHHmmss');

fclose(fid);
save(full, 'C');