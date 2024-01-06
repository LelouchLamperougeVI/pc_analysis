%% supplementary fig1 - lesion vs control from Esteves et al. 2021 jneuro
clear all

% make them damn objects
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
    obj(ii) = ensemble(abfs{ii}, 'planes', 2);
%     obj(ii).rm_redund;
    obj(ii).perform_analysis;
    obj(ii).bayes_infer2;
end

save('/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation/obj.mat', 'obj');


%%
clear all
load('/mnt/storage/esteves_et_al_jneuro2021_for_supp1_in_chang_et_al_m2_reactivation/obj.mat');

% partition analysis structs
l_analysis = []; l_err = []; l_cm = [];
c_analysis = []; c_err = []; c_cm = [];
for ii = 1:length(obj)
    temp = regexpi(obj(ii).session.animal, '^hl*');
    cm = [discretize(obj(ii).bayes.x, obj(ii).bayes.xbins), discretize(obj(ii).bayes.decoded, obj(ii).bayes.xbins)];
    cm(any(isnan(cm), 2), :) = [];
    cm = accumarray(cm, 1, [obj(ii).bayes.ops.bins, obj(ii).bayes.ops.bins]);
    if isempty(temp)
        c_analysis = [c_analysis; obj(ii).analysis];
        c_err = [c_err, obj(ii).bayes.err];
        c_cm = cat(3, c_cm, cm);
    else
        l_analysis = [l_analysis; obj(ii).analysis];
        l_err = [l_err, obj(ii).bayes.err];
        l_cm = cat(3, l_cm, cm);
    end
end

c_cm = sum(c_cm, 3); c_cm = c_cm ./ sum(c_cm, 2); % P(D | x)
l_cm = sum(l_cm, 3); l_cm = l_cm ./ sum(l_cm, 2);

% panel a - stacks
c_stack = arrayfun(@(x) x.stack(:, x.pc_list), c_analysis, 'UniformOutput', false);
l_stack = arrayfun(@(x) x.stack(:, x.pc_list), l_analysis, 'UniformOutput', false);
c_stack = cat(2, c_stack{:});
l_stack = cat(2, l_stack{:});

figure
subplot(1, 2, 1);
[~, idx] = max(c_stack, [], 1);
[~, idx] = sort(idx);
imagesc(-c_stack(:, idx)')
rbmap('caxis', [-1 0], 'interp', 255);

subplot(1, 2, 2);
[~, idx] = max(l_stack, [], 1);
[~, idx] = sort(idx);
imagesc(l_stack(:, idx)')
rbmap('caxis', [0 1], 'interp', 255);

% panel b - place field location densities
c_loc = arrayfun(@(x) cat(1, x.width{x.pc_list}), c_analysis, 'UniformOutput', false);
c_loc = cat(1, c_loc{:});
c_loc = c_loc(:, 2);
l_loc = arrayfun(@(x) cat(1, x.width{x.pc_list}), l_analysis, 'UniformOutput', false);
l_loc = cat(1, l_loc{:});
l_loc = l_loc(:, 2);

nboot = 1e3;
edges = 0.5:50.5;
ci = zeros(length(edges) - 1, 2, nboot);
for ii = 1:nboot
    temp = randsample(c_loc, length(c_loc), true);
    ci(:, 1, ii) = histcounts(temp, edges, 'Normalization', 'probability');
    temp = randsample(l_loc, length(l_loc), true);
    ci(:, 2, ii) = histcounts(temp, edges, 'Normalization', 'probability');
end
ci = prctile(ci, [2.5, 97.5], 3);
loc = cat(2, histcounts(c_loc, edges, 'Normalization', 'probability')', histcounts(l_loc, edges, 'Normalization', 'probability')');
g = repmat((1:2), [size(loc, 1), 1]);
x = repmat(conv(edges, [.5 .5], "valid")', [1 2]);

g = gramm('x', x(:), 'y', loc(:), 'color', g(:), 'ymin', reshape(ci(:, :, 1), [], 1), 'ymax', reshape(ci(:, :, 2), [], 1));
g.geom_line();
g.geom_interval();
g.axe_property('ylim', [0 .1]);
figure
g.draw();

% panel c - bayesian decoding
x = linspace(0, 150, size(l_err, 1));
x = repmat(x, 1, 2);
y = cat(2, mean(c_err'), mean(l_err'));
ci = cat(2, sem(c_err'), sem(l_err'));
g = repelem((1:2), size(l_err, 1));
g = gramm('x', x, 'y', y, 'ymin', y - ci, 'ymax', y + ci, 'color', g);
g.geom_line();
g.geom_interval();
g.axe_property('ylim', [0 35]);
figure
g.draw();

figure
subplot(1, 2, 1)
imagesc(imgaussfilt(c_cm, .5))
axis square
colormap plasma
colorbar
caxis([0 .6])
xlabel('decoded'); ylabel('real'); title('P(D | x)');
subplot(1, 2, 2)
imagesc(imgaussfilt(l_cm, .5))
axis square
colormap plasma
caxis([0 .6])
colorbar

% panel d - number of place fields

c_num_pf = arrayfun(@(x) cellfun(@(y) size(y, 1), x.width(x.pc_list)), c_analysis, 'UniformOutput', false);
c_num_pf = cat(2, c_num_pf{:});
l_num_pf = arrayfun(@(x) cellfun(@(y) size(y, 1), x.width(x.pc_list)), l_analysis, 'UniformOutput', false);
l_num_pf = cat(2, l_num_pf{:});

g = repelem((1:2)', [length(c_num_pf); length(l_num_pf)]);
num_pf = cat(2, c_num_pf, l_num_pf);

g = gramm('x', num_pf, 'color', g);
g.stat_bin('normalization', 'probability', 'edges', .5:1:3.5);
g.axe_property('ylim', [0 1]);
figure
g.draw();

p = chisq2(c_num_pf, l_num_pf)


%% table for one-cue/off-cue decoding error
c_cue_err = [mean(c_err(belt_idx, :)); mean(c_err(~belt_idx, :))]';
l_cue_err = [mean(l_err(belt_idx, :)); mean(l_err(~belt_idx, :))]';
cue_err = cat(1, c_cue_err, l_cue_err);
cue_err = arrayfun(@(x) x, cue_err, 'UniformOutput', false);
g = repelem({'control'; 'lesion'}, [size(c_cue_err, 1), size(l_cue_err, 1)], 2);
cue = repmat({'on', 'off'}, size(cue_err, 1), 1);
subject = repmat((1:size(cue_err, 1))', 1, 2);
subject = arrayfun(@(x) x, subject, 'UniformOutput', false);
tab = [subject(:), g(:), cue(:), cue_err(:)];
tab = cell2table(tab, 'VariableNames', {'subject', 'sx', 'cue', 'error'});

% writetable(tab, '/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/R/cue_decoding_error.csv');
