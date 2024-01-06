%% Make them damn objects
clear all

root = '/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation';
animal = dir(root);

abfs = {};
for a = 3:length(animal)
    date = dir(fullfile(root, animal(a).name));
    for d = 3:length(date)
        session = dir(fullfile(root, animal(a).name, date(d).name, '*.abf'));
        abfs = cat(2, abfs, fullfile(root, animal(a).name, date(d).name, {session.name}));
    end
end
abfs = abfs';

for ii = 1:length(abfs)
    obj(ii) = ensemble(abfs{ii}, 'planes', 1);
%     obj(ii).rm_redund;
    obj(ii).perform_analysis;
    obj(ii).bayes_infer2;
end

% save('/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation/obj.mat', 'obj');


%% Index control/lesion
lesioned = false(length(obj), 1);
for ii = 1:length(obj)
    temp = regexpi(obj(ii).session.animal, '^hl*');
    lesioned(ii) = ~isempty(temp);
end

%% Make error csv for decoding
err = arrayfun(@(x) x.bayes.err, obj, 'UniformOutput', false);
err = cell2mat(err);

err_tab = array2table([lesioned, err'], 'VariableNames', cat(2, {'lesion'}, arrayfun(@(x) ['x', num2str(x)], 1:size(err, 1), 'UniformOutput', false)));
% writetable(err_tab, '/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/R/err.csv');