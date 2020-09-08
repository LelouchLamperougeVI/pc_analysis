clear all

root = '/mnt/storage/rrr_magnum/M2';
animal = 'EE003';
date = '2018_12_29';
session = '3';

ass = ensemble(fullfile(root, animal, date, session, [date '_' session '.abf']));

%%
ass.set_ops('e_size',10);
ass.set_ops('clust_method','thres');
ass.cluster;
ass.set_ops('order','pc');

%% plot rest-run
figure
subplot(1, 2, 1);
deconv = ass.analysis.deconv(:, ass.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(-deconv');
colormap gray
subplot(1, 2, 2);
deconv = ass.twop.deconv(:, ass.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, ass.ops.sig * ass.twop.fs);
imagesc(-deconv');
colormap gray


%% Data tree
clear all

date_fcn = @(x) datetime(x, 'inputformat', 'yyyy_MM_dd');
root = '/mnt/storage/rrr_magnum/M2';
% tree = dirTree(root, 'isdir', [1 1 1 1 1 0], 're', {nan, nan, nan, '3', 'plane.', 'deconv.mat'}, 'fcn', {nan,nan,date_fcn,nan,nan,nan});
tree = dirTree(root, 'isdir', [1 1 1 0], 're', {nan, nan, nan, '_2.abf'}, 'fcn', {nan,nan,date_fcn,nan});
anal_list = tree.fullfile;
% anal_list( cellfun(@length, anal_list) < 60 ) = [];

exceptions = {'/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_2.abf';
    '/mnt/storage/rrr_magnum/M2/EE001/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC003/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_20/2018_12_20_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_21/2018_12_21_2.abf';
    '/mnt/storage/rrr_magnum/M2/EC002/2018_12_28/2018_12_28_2.abf'};

anal_list = setdiff(anal_list, exceptions);
clearvars -except anal_list

%% make batch analysis
% for ii = 8:length(anal_list)
for ii = 54:length(anal_list)
    lfp = LFP(anal_list{ii});
    lfp.perform_analysis;
    lfp.save('analysis');
    lfp.save('chan');
end


%% check analysis
clear all

root = '/mnt/storage/rrr_magnum/M2';
tree = dirTree(root, 'isdir', [1 1 1 0], 're', {nan, nan, nan, 'analysis.mat'}, 'fcn', {nan,nan,nan,nan});
anal_list = tree.fullfile;
anal_list( cellfun(@length, anal_list) < 56 ) = [];

for ii = 1:length(anal_list)
    load(anal_list{ii});
    disp([anal_list{ii} ': ' num2str(length(analysis.pc_list) / size(analysis.original_deconv, 2))]);
end


%% Aubrey's animals with mutable cues
stack1 = run1.analysis.stack;
stack2 = run2.analysis.stack;

[~, order] = max(stack1, [], 1);
[~, order] = sort(order);
order = intersect(order, run1.analysis.pc_list, 'stable');

figure
subplot(2, 2, 1);
imagesc(stack1(:, order)')
colormap jet
subplot(2, 2, 2);
imagesc(stack2(:, order)')
colormap jet

[~, order] = max(stack2, [], 1);
[~, order] = sort(order);
order = intersect(order, run2.analysis.pc_list, 'stable');

subplot(2, 2, 3);
imagesc(stack1(:, order)')
colormap jet
subplot(2, 2, 4);
imagesc(stack2(:, order)')
colormap jet


%% Fig1C
cd /home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB

run = LFP('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_2.abf');
run.perform_analysis;

rest1 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_1.abf');
rest2 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');

rest1.import_analysis(run.analysis)
rest1.set_ops('order','pc');
rest2.import_analysis(run.analysis)
rest2.set_ops('order','pc');

%%
figure
ax1(1) = subplot(5, 3, [1 4 7]);
deconv = rest1.twop.deconv(:, rest1.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest1.ops.sig * rest1.twop.fs);
imagesc('xdata', rest1.twop.ts, 'cdata', -deconv');
colormap gray
ax1(2) = subplot(5, 3, 10);
plot(rest1.behavior.ts, rest1.behavior.speed);
ylim([-20 60])
ax1(3) = subplot(5, 3, 13);
plot(rest1.behavior.ts, rest1.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax1, 'x')
xlim([300 400])

ax2(1) = subplot(5, 3, [2 5 8]);
deconv = run.twop.deconv(:, run.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest1.ops.sig * run.twop.fs);
imagesc('xdata', run.twop.ts, 'cdata', -deconv');
colormap gray
ax2(2) = subplot(5, 3, 11);
plot(run.behavior.ts, run.behavior.speed);
ylim([-20 60])
ax2(3) = subplot(5, 3, 14);
plot(run.behavior.ts, run.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax2, 'x')
xlim([200 300])

ax3(1) = subplot(5, 3, [3 6 9]);
deconv = rest2.twop.deconv(:, rest2.analysis.order);
deconv = (deconv - min(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, rest2.ops.sig * rest2.twop.fs);
imagesc('xdata', rest2.twop.ts, 'cdata', -deconv');
colormap gray
ax3(2) = subplot(5, 3, 12);
plot(rest2.behavior.ts, rest2.behavior.speed);
ylim([-20 60])
ax3(3) = subplot(5, 3, 15);
plot(rest2.behavior.ts, rest2.behavior.pos_raw);
ylim([-10 160])
linkaxes(ax3, 'x')
xlim([300 400])

%% Fig1D
mimg = rest2.topo.mimg;
mimg = (mimg - min(mimg(:))) ./ range(mimg(:));
mimg = mimg ./ .25; %adjust contrast
mimg(mimg>1) = 1;
mask = rest2.topo.maskNeurons;

f = rest2.twop.deconv(rest2.twop.ts > 310.5 & rest2.twop.ts < 311.5, :);
f = mean(f) ./ range(rest2.twop.deconv);
im = repmat(mimg, 1, 1, 3);
c = rbmap('caxis',[0 max(f)]);
f = discretize(f, size(c, 1));
for ii = 1:length(f)
    [x, y] = find(mask == ii);
    for jj = 1:length(x)
        im(x(jj), y(jj), :) = permute(c(f(ii), :), [1 3 2]);
    end
end

imwrite(im, 'fig1d_mask_react.png')
% figure
% imshow(im);

[~, idx] = sort(f, 'descend');
[~, idx] = max( run.analysis.stack(:, idx(1:5)), [], 1 );
idx = round(median(idx));

f = run.analysis.stack(idx-1:idx+1, :);
f = mean(f);
im = repmat(mimg, 1, 1, 3);
c = rbmap('caxis',[0 max(f)]);
f = discretize(f, size(c, 1));
for ii = 1:length(f)
    [x, y] = find(mask == ii);
    for jj = 1:length(x)
        im(x(jj), y(jj), :) = permute(c(f(ii), :), [1 3 2]);
    end
end

imwrite(im, 'fig1d_mask_run.png')
% figure
% imshow(im);

%% Fig2
rest2.set_ops('e_size',5);
rest2.set_ops('clust_method','thres');
rest2.set_ops('sig', .2);
rest2.remove_mvt;
rest2.cluster;
rest2.topography;
rest2.detect_sce;

% rest2.plot('clust_topo')
rest2.set_ops('order','cluster')

rest2.plot('sce')
rest2.plot('clust_corr')
rest2.plot('tree')
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{4})
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{5})
plot_analysis(run.analysis, [1 0 0 ], rest2.ensembles.clust{7})

plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{4})
plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{5})
plot_analysis(run.analysis, [0 1 0 ], rest2.ensembles.clust{7})
% nice clusts: 4, 5, 7

%%
figure
deconv = fast_smooth(rest2.twop.deconv, rest2.ops.sig * rest2.twop.fs);
deconv = (deconv - min(deconv, [], 'omitnan')) ./ range(deconv);
imagesc('xdata', rest2.twop.ts, 'cdata', -deconv(:, rest2.ensembles.clust_order)');
colormap gray
xlim([355 415])



%% Batch make/save analysis old M2
chan = [1; 2; 3; 5; NaN; NaN; NaN];

root = '/mnt/storage/HaoRan/RRR_motor/M2';

animals = dir(fullfile(root, 'RSC*'));
animals = {animals.name};

for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        disp(['current session: ' fullfile(animals{a}, sessions{s})])
        save(fullfile(root, animals{a}, sessions{s}, 'channels.mat'), 'chan');
        lfp = LFP(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_2.abf']));
        lfp.perform_analysis;
        lfp.save('analysis');
    end
end

%% Batch make/save analysis new M2
clear all

root = '/mnt/storage/rrr_magnum/M2';

animals = dir(fullfile(root, 'E*'));
animals = {animals.name};

for a = 4:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        lfp = LFP(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_2.abf']));
        lfp.perform_analysis;
        lfp.save('analysis');
    end
end


%% Rest analysis
clear all

traj_thres = .2;

% root = '/mnt/storage/HaoRan/RRR_motor/M2';
root = '/mnt/storage/rrr_magnum/M2';

% animals = dir(fullfile(root, 'RSC*'));
animals = dir(fullfile(root, 'E*'));
animals = {animals.name};

EV = [];
frac = []; %[rest1 rest2 run]
frac_clust = cell(2,1);
clust_stacks = cell(2,1);
trajectories = cell(2,1);
loc = {};
loc_clust = cell(2,1);
swr_stack = cell(2,1);
clusts = cell(2,1);
SI = [];
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.set_ops('e_size',5);
        rest1.set_ops('clust_method','thres');
        rest1.set_ops('sig', .2);
        rest1.remove_mvt;
        rest1.cluster;
        rest1.swr_window;
%         rest1.topography;
%         rest1.detect_sce;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
        rest2.swr_window;
%         rest2.topography;
%         rest2.detect_sce;
        [temp1, temp2] = ev(rest1, rest1.analysis.original_deconv, rest2);
        EV = cat(1, EV, [temp1 temp2]);
        
        frac = cat(1, frac, [length(intersect( cell2mat(rest1.ensembles.clust), rest1.analysis.pc_list )) / length(cell2mat(rest1.ensembles.clust)), ...
        length(intersect( cell2mat(rest2.ensembles.clust), rest1.analysis.pc_list )) / length(cell2mat(rest2.ensembles.clust)), ...
        length(rest1.analysis.pc_list) / size(rest1.twop.deconv, 2)]);
    
        swr_stack{1} = cat(1, swr_stack{1}, {rest1.ensembles.swr.all});
        swr_stack{2} = cat(1, swr_stack{2}, {rest2.ensembles.swr.all});
        clusts{1} = cat(2, clusts{1}, {rest1.ensembles.clust});
        clusts{2} = cat(2, clusts{2}, {rest2.ensembles.clust});
        SI = cat(2, SI, {rest1.analysis.SI});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        for c = 1:length(rest1.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest1.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            frac_clust{1} = cat(1, frac_clust{1}, [length(list) length(list)/length(rest1.ensembles.clust{c})]);
%             clust_stacks{1} = cat(2, clust_stacks{1}, mean(stack, 2));
            clust_stacks{1} = cat(1, clust_stacks{1}, {stack});
            trajectories{1} = cat(2, trajectories{1}, any(stack > traj_thres, 2));
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({rest1.analysis.width{~cellfun(@isempty, rest1.analysis.width)}}');
        loc{end+1} =  temp(:, 2);
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{1} = cat(1, loc_clust{1}, {temp});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        for c = 1:length(rest2.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest2.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            frac_clust{2} = cat(1, frac_clust{2}, [length(list) length(list)/length(rest2.ensembles.clust{c})]);
%             clust_stacks{2} = cat(2, clust_stacks{2}, mean(stack, 2, 'omitnan'));
            clust_stacks{2} = cat(1, clust_stacks{2}, {stack});
            trajectories{2} = cat(2, trajectories{2}, any(stack > traj_thres, 2));
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{2} = cat(1, loc_clust{2}, {temp});
    end
end


%% Combine run for EE003/2018_12_19
lfp1 = LFP('/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_2.abf');
behav1 = lfp1.behavior;
lfp2 = LFP('/mnt/storage/rrr_magnum/M2/EE003/2018_12_19/2018_12_19_3.abf');
behav2 = lfp2.behavior;
ts1 = lfp1.twop.ts;
ts2 = lfp2.twop.ts;
dec1 = lfp1.twop.deconv;
dec2 = lfp2.twop.deconv;
ts = cat(1, ts1, ts2 + ts1(end) + median(diff(ts1)));
dec = cat(1, dec1, dec2);

beh = behav1;
beh.speed_raw_noSmooth = [beh.speed_raw_noSmooth behav2.speed_raw_noSmooth];
beh.speed_raw = [beh.speed_raw; behav2.speed_raw];

beh.ts = [beh.ts behav2.ts + ts1(end) + median(diff(ts1))];

beh.pos_norm = [beh.pos_norm behav2.pos_norm];
beh.trial= [beh.trial behav2.trial+max(beh.trial)];
beh.speed = [beh.speed behav2.speed];
beh.pos_cum = [beh.pos_cum behav2.pos_cum];
beh.pos_raw = [beh.pos_raw behav2.pos_raw];

tcs.tt = ts';
[beh, dec] = convert_behavior(beh, tcs, dec);
analysis = pc_batch_analysis(beh, dec);


%%
clear all

frac_thres = 3;

frac_clust1 = load('/mnt/storage/HaoRan/RRR_motor/M2/frac_clust.mat');
frac_clust1 = frac_clust1.frac_clust;
frac_clust2 = load('/mnt/storage/rrr_magnum/M2/frac_clust.mat');
frac_clust2 = frac_clust2.frac_clust;

trajectories1 = load('/mnt/storage/HaoRan/RRR_motor/M2/traj.mat');
trajectories1 = trajectories1.trajectories;
trajectories2 = load('/mnt/storage/rrr_magnum/M2/traj.mat');
trajectories2 = trajectories2.trajectories;

trajectories = [trajectories1{1} trajectories2{1}];
frac_clust = [frac_clust1{1}(:,1); frac_clust2{1}(:,1)];
trajectories = trajectories(:, frac_clust >= frac_thres);
[l1, s1, e1] = traj_length(trajectories);
% l = cellfun(@(x) any(x > 20), l);

trajectories = [trajectories1{2} trajectories2{2}];
frac_clust = [frac_clust1{2}(:,1); frac_clust2{2}(:,1)];
trajectories = trajectories(:, frac_clust >= frac_thres);
[l2, s2, e2] = traj_length(trajectories);
% l = cellfun(@(x) any(x > 20), l);

belt;

figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

figure
histogram(cell2mat(s1), 50)
hold on
histogram(cell2mat(e1), 50)
plot(linspace(0, 150, 50), belt_idx.*10)

figure
histogram(cell2mat(s2), 50)
hold on
histogram(cell2mat(e2), 50)
plot(linspace(0, 150, 50), belt_idx.*10)


%% Figure 3
clear all

circ_dist = @(X0, X) min(cat(3, mod(X0 - X, 50), mod(X - X0, 50)), [], 3) .* 150 ./ 50;

thres = 3; %min number of pc per ensemble

loc1 = load('/mnt/storage/HaoRan/RRR_motor/M2/loc.mat');
loc2 = load('/mnt/storage/rrr_magnum/M2/loc.mat');
loc = cat(2, loc1.loc, loc2.loc);
loc_clust{1} = cat(1, loc1.loc_clust{1}, loc2.loc_clust{1});
loc_clust{2} = cat(1, loc1.loc_clust{2}, loc2.loc_clust{2});

d1 = []; d2 = [];
loc_cum1 = []; loc_cum2 = [];

s1 = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{1}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = cellfun(@(x) accumarray(x + 1, 1), loc_clust{1}, 'uniformoutput', false);
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 1:length(g{ii})
        d = loc{ii}(loc_clust{1}{ii} == g{ii}(jj));
        loc_cum1 = cat(1, loc_cum1, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d1 = cat(1, d1, mean(d));
    end
end
g = arrayfun(@(x) ismember(loc_clust{1}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s1 = arrayfun(@(x) s1{x}(g{x}), 1:length(s1), 'uniformoutput', false);

s2 = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, circ_dist), 1:length(loc), 'uniformoutput', false);
% s = arrayfun(@(x) silhouette(loc{x}, loc_clust{2}{x}, 'Euclidean'), 1:length(loc), 'uniformoutput', false);
g = cellfun(@(x) accumarray(x + 1, 1), loc_clust{2}, 'uniformoutput', false);
g = cellfun(@(x) find(x(2:end) >= thres), g, 'uniformoutput', false);
for ii = 1:length(g)
    for jj = 1:length(g{ii})
        d = loc{ii}(loc_clust{2}{ii} == g{ii}(jj));
        loc_cum2 = cat(1, loc_cum2, d);
        d = circ_dist(d, d');
        d = d(logical(triu(ones(size(d)), 1)));
        d2 = cat(1, d2, mean(d));
    end
end
g = arrayfun(@(x) ismember(loc_clust{2}{x}, g{x}), 1:length(loc), 'uniformoutput', false);
s2 = arrayfun(@(x) s2{x}(g{x}), 1:length(s2), 'uniformoutput', false);

figure
cdfplot(cell2mat(s1'));
hold on
cdfplot(cell2mat(s2'));

[~, p] = kstest2(cell2mat(s1'), cell2mat(s2'));
disp(['kstest silhouette: ' num2str(p)])

figure
boxplot([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)])

p = kruskalwallis([cell2mat(s1'); cell2mat(s2')], [ones(length(cell2mat(s1')), 1); 2 .* ones(length(cell2mat(s2')), 1)]);
disp(['kruskal-wallis silhouette: ' num2str(p)])

figure;
violin({d1, d2}, 'labels', {'rest1', 'rest2'})
% violin({d1, d2}, 'labels', {'rest1', 'rest2'}, 'scatter')

p = kruskalwallis([d1; d2], [ones(length(d1), 1); 2 .* ones(length(d2), 1)]);
disp(['kruskal-wallis dist: ' num2str(p)])

figure;
histogram(loc_cum1, 50, 'normalization', 'probability')
hold on
histogram(loc_cum2, 50, 'normalization', 'probability')


%% Trajectories
clear all

clust_stacks1 = load('/mnt/storage/HaoRan/RRR_motor/M2/clust_stacks.mat');
clust_stacks2 = load('/mnt/storage/rrr_magnum/M2/clust_stacks.mat');
clust_stacks{1} = cat(1, clust_stacks1.clust_stacks{1}, clust_stacks2.clust_stacks{1});
clust_stacks{2} = cat(1, clust_stacks1.clust_stacks{2}, clust_stacks2.clust_stacks{2});

fr_thres = .5;
traj_thres = 3; %min number of pc per ensemble

l1 = cell(length(clust_stacks{1}), 1); s1=l1; e1=l1;
for c = 1:length(clust_stacks{1})
    stack = clust_stacks{1}{c};
    traj = any(stack > fr_thres, 2);
    [~, starts, ends] = traj_length(traj, 1);

    stack = repmat(stack, 2, 1);
    idx = false(length(starts{1}), 1);
    for t = 1:length(starts{1})
        temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
        temp = any(temp > fr_thres, 1);
        idx(t) = sum(temp) < traj_thres;
    end

    [l1{c}, s1{c}, e1{c}] = traj_length(traj);
    l1{c} = l1{c}{1}(~idx);
    s1{c} = s1{c}{1}(~idx);
    e1{c} = e1{c}{1}(~idx);
end

l2 = cell(length(clust_stacks{2}), 1); s2=l2; e2=l2;
for c = 1:length(clust_stacks{2})
    stack = clust_stacks{2}{c};
    traj = any(stack > fr_thres, 2);
    [~, starts, ends] = traj_length(traj, 1);

    stack = repmat(stack, 2, 1);
    idx = false(length(starts{1}), 1);
    for t = 1:length(starts{1})
        temp = stack(starts{1}(t) : (ends{1}(t) - 1 + length(traj) * (starts{1}(t) > (ends{1}(t)))), :);
        temp = any(temp > fr_thres, 1);
        idx(t) = sum(temp) < traj_thres;
    end

    [l2{c}, s2{c}, e2{c}] = traj_length(traj);
    l2{c} = l2{c}{1}(~idx);
    s2{c} = s2{c}{1}(~idx);
    e2{c} = e2{c}{1}(~idx);
end

belt;

figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

figure
histogram(cell2mat(s1), 50, 'normalization', 'probability')
title('rest1 start')
xlim([0 150]); ylim([0 .2])
figure
histogram(cell2mat(e1), 50, 'normalization', 'probability')
title('rest1 end')
xlim([0 150]); ylim([0 .2])
% plot(linspace(0, 150, 50), belt_idx.*.01)

figure
histogram(cell2mat(s2), 50, 'normalization', 'probability')
title('rest2 start')
xlim([0 150]); ylim([0 .2])
figure
histogram(cell2mat(e2), 50, 'normalization', 'probability')
title('rest2 end')
xlim([0 150]); ylim([0 .2])
% plot(linspace(0, 150, 50), belt_idx.*.01)


bins = 25;

edges = linspace(0, 150, bins+1);
P_se = accumarray([discretize(cell2mat(s1), edges) discretize(cell2mat(e1), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_s_cond_e = P_se ./ sum(P_se, 1);
P_e_cond_s = P_se ./ sum(P_se, 2);

% figure
% imagesc(P_se)
% colormap jet
% colorbar
% rbmap('caxis', [0 1]);
% axis image
figure
imagesc(P_e_cond_s)
rbmap('caxis', [0 1]);
colorbar
caxis([0 1])
axis image
% figure
% imagesc(P_s_cond_e)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image


P_se = accumarray([discretize(cell2mat(s2), edges) discretize(cell2mat(e2), edges)], 1, [bins bins]);
P_se = P_se ./ sum(P_se(:));
P_s_cond_e = P_se ./ sum(P_se, 1);
P_e_cond_s = P_se ./ sum(P_se, 2);

% figure
% imagesc(P_se)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image
figure
imagesc(P_e_cond_s)
rbmap('caxis', [0 1]);
colorbar
caxis([0 1])
axis image
% figure
% imagesc(P_s_cond_e)
% rbmap('caxis', [0 1]);
% colorbar
% caxis([0 1])
% axis image


%% Fig3c
figure
% plot(sum(trajectories'))
cdfplot(cell2mat(l1));
hold on
cdfplot(cell2mat(l2))

n = 464;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 465;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 353;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 334;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 26;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 68;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone
n = 200;
figure; [~, idx] = max(clust_stacks{2}{n}); [~, idx] = sort(idx); imagesc(-clust_stacks{2}{n}(:, idx)'); title([num2str(n) '; length:' num2str(l2{n}) '; start:' num2str(s2{n}) '; end:' num2str(e2{n})]); colormap bone


%% Reconstructed trajectories
P_e_cond_s(isnan(P_e_cond_s)) = 0;
P_se(isnan(P_se)) = 0;

% g = digraph(P_e_cond_s);
g = digraph(P_se);
LWidths = 10*g.Edges.Weight/max(g.Edges.Weight);
LWidths(isnan(LWidths)) = 0;

figure
plot(g, 'layout', 'circle', 'LineWidth', LWidths, 'edgecdata', g.Edges.Weight, 'arrowsize', 0)
axis square
% rbmap('caxis',[0 1]);


%% bad LFPs
bad: EC002, EE003
really good: EC007, EE006


%%
clear all
load('/mnt/storage/rrr_magnum/M2/swr_stack.mat');

sig = 1; % smoothing factor

n_stack = []; % ensemble neurons stack
e_stack = []; % ensemble stack
p_stack = []; % population stack

% for s = 1:length(clusts{1})
%     if size(swr_stack{1}{s}, 1) ~= 39; continue; end
%     for e = 1:length(clusts{1}{s})
%         temp = swr_stack{1}{s}(:, clusts{1}{s}{e}, :);
%         temp = mean(temp, 3, 'omitnan');
%         n_stack = cat(2, n_stack, temp);
%         temp = mean(temp, 2);
%         e_stack = cat(2, e_stack, temp);
%     end
%     temp = swr_stack{1}{s};
%     temp = mean(temp, 3, 'omitnan');
%     temp = mean(temp, 2);
%     p_stack = cat(2, p_stack, temp);
% end

for s = 1:length(clusts{2})
    if size(swr_stack{2}{s}, 1) ~= 39; continue; end
    for e = 1:length(clusts{2}{s})
        temp = swr_stack{2}{s}(:, clusts{2}{s}{e}, :);
        temp = mean(temp, 3, 'omitnan');
        n_stack = cat(2, n_stack, temp);
        temp = mean(temp, 2);
        e_stack = cat(2, e_stack, temp);
    end
    temp = swr_stack{2}{s};
    temp = mean(temp, 3, 'omitnan');
    temp = mean(temp, 2);
    p_stack = cat(2, p_stack, temp);
end

% n_stack = n_stack(10:29,:);
% e_stack = e_stack(10:29,:);
% p_stack = p_stack(10:29,:);

[~, idx] = max(e_stack);
[~, idx] = sort(idx);
figure;
temp = e_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('e_stack')
caxis([-.1 .1])

[~, idx] = max(n_stack);
[~, idx] = sort(idx);
figure;
temp = n_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('n_stack')
caxis([-.1 .1])

[~, idx] = max(p_stack);
[~, idx] = sort(idx);
figure;
temp = p_stack(:,idx);
temp = fast_smooth(temp, sig);
imagesc(temp');
rbmap('caxis',[-.1 .1]);
set(gca, 'PlotBoxAspectRatio', [.5 1 1]);
title('p_stack')
caxis([-.1 .1])

% temp = fast_smooth(e_stack, sig);
% h = errorshade(mean(temp,2,'omitnan'), sem(temp, 2));
% 
% temp = fast_smooth(p_stack, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2), 'h', h)
% 
% temp = cellfun(@(x) mean(x, 3, 'omitnan'), swr_stack{1}, 'uniformoutput', false);
% temp = cell2mat(temp(setdiff(1:68, 35))');
% temp = fast_smooth(temp, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2))
% 
% temp = cellfun(@(x) mean(x, 3, 'omitnan'), swr_stack{2}, 'uniformoutput', false);
% temp = cell2mat(temp(setdiff(1:68, 35))');
% temp = fast_smooth(temp, sig);
% errorshade(mean(temp,2,'omitnan'), sem(temp, 2))


%% Fig5a
clear all

animal = 'EC007';
date = '2019_05_29';

lfp = ensemble(fullfile('/mnt/storage/rrr_magnum/M2/', animal, date, [date '_3.abf']));
lfp.set_ops('e_size',5);
lfp.set_ops('clust_method','thres');
lfp.set_ops('sig', .2);
lfp.remove_mvt;
lfp.cluster;
lfp.set_ops('order','cluster')

% lfp.swr_window;
% lfp.detect_sce;
% lfp.sce_spectrum;

deconv = lfp.twop.deconv(:, lfp.ensembles.clust_order);
deconv = (deconv - nanmin(deconv)) ./ range(deconv);
deconv = fast_smooth(deconv, lfp.ops.sig * lfp.twop.fs);

figure
ax(1) = subplot(6,1,1:2);
imagesc('xdata', lfp.twop.ts, 'cdata', -deconv')
colormap gray
ax(2) = subplot(6,1,3);
plot(lfp.lfp.ts, lfp.lfp.lfp);
hold on
plot(lfp.lfp.swr.swr_on, max(lfp.lfp.lfp) .* ones(length(lfp.lfp.swr.swr_on), 1), 'k*');
ax(3) = subplot(6,1,4);
plot(lfp.lfp.ts, lfp.lfp.swr.swr);
ax(4) = subplot(6,1,5);
plot(lfp.lfp.ts, lfp.lfp.gamma);
ax(5) = subplot(6,1,6);
plot(lfp.behavior.ts, lfp.behavior.speed);
linkaxes(ax, 'x')














