clear all

root = '/mnt/storage/HaoRan/RRR_motor/M2';
animals = dir(fullfile(root, 'RSC*'));

% root = '/mnt/storage/rrr_magnum/M2';
% animals = dir(fullfile(root, 'E*'));

animals = {animals.name};

count = 1;

f_lin = linspace(0, 10, 2^8);
f = logspace(-2, 1, 2^8);
p = zeros(length(f), 3);
p_lin = p;
p_neur = p;
pfr = zeros(3, 1);
R = [];
coh = zeros(length(f), 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    disp(['running ' num2str(a) '/' num2str(length(animals))]);
    
    for s = 1:length(sessions)
        load(fullfile(root, animals{a}, sessions{s}, 'analysis.mat'));
        vel = analysis.behavior.unit_vel;
        vel(isnan(vel)) = 0;
        for t = 1:3
            path = fullfile(root, animals{a}, sessions{s}, num2str(t));
            flag = dir(fullfile(path, 'plane1'));
            if isempty(flag)
                path = fullfile(path, 'Plane1');
            else
                path = fullfile(path, 'plane1');
            end
            load(fullfile(path, 'deconv.mat'));
            load(fullfile(path, 'timecourses.mat'));
            
            fs = 1 / median(diff(tcs.tt));
            d = ca_filt(deconv);
            d = mean(d > 0, 2);
            d = d .* fs; % normalise to mean #spikes per second
            %             d = mean(deconv, 2);
            pfr(t, count) = mean(d);
            
            wdw = round(fs*200);
            nover = round(wdw * .75);
            p(:, t, count) = pwelch(d, hann(wdw), nover, f, fs, 'psd');
            p_lin(:, t, count) = pwelch(d, hann(wdw), nover, f_lin, fs, 'psd');
            
            neur = mean((tcs.neuropil - tcs.baseline) ./ tcs.baseline .* 100, 2);
            p_neur(:, t, count) = pwelch(neur, hann(wdw), nover, f, fs, 'psd');
            
            if t == 2
                coh(:, count) = sqrt(mscohere(d(:), vel(:), hann(wdw), nover, f, fs));
            end
            
            deconv=fast_smooth(deconv, .2 * fs);
            deconv=(deconv-mean(deconv))./std(deconv); %zscore
            R{count}(:, :, t) = corr(deconv);
        end
        disp([num2str(s) '/' num2str(length(sessions))]);
        count = count + 1;
    end
end


root = '/mnt/storage/rrr_magnum/M2';
animals = dir(fullfile(root, 'E*'));

animals = {animals.name};


for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    disp(['running ' num2str(a) '/' num2str(length(animals))]);
    
    for s = 1:length(sessions)
        load(fullfile(root, animals{a}, sessions{s}, 'analysis.mat'));
        vel = analysis.behavior.unit_vel;
        vel(isnan(vel)) = 0;
        for t = 1:3
            path = fullfile(root, animals{a}, sessions{s}, num2str(t));
            flag = dir(fullfile(path, 'plane1'));
            if isempty(flag)
                path = fullfile(path, 'Plane1');
            else
                path = fullfile(path, 'plane1');
            end
            load(fullfile(path, 'deconv.mat'));
            load(fullfile(path, 'timecourses.mat'));
            
            fs = 1 / median(diff(tcs.tt));
            d = ca_filt(deconv);
            d = mean(d > 0, 2);
            d = d .* fs; % normalise to mean #spikes per second
            %             d = mean(deconv, 2);
            pfr(t, count) = mean(d);
            
            wdw = round(fs*200);
            nover = round(wdw * .75);
            p(:, t, count) = pwelch(d, hann(wdw), nover, f, fs, 'psd');
            p_lin(:, t, count) = pwelch(d, hann(wdw), nover, f_lin, fs, 'psd');
            
            neur = mean((tcs.neuropil - tcs.baseline) ./ tcs.baseline .* 100, 2);
            p_neur(:, t, count) = pwelch(neur, hann(wdw), nover, f, fs, 'psd');
            
            if t == 2
                coh(:, count) = sqrt(mscohere(d(:), vel(:), hann(wdw), nover, f, fs));
            end
            
            deconv=fast_smooth(deconv, .2 * fs);
            deconv=(deconv-mean(deconv))./std(deconv); %zscore
            R{count}(:, :, t) = corr(deconv);
        end
        disp([num2str(s) '/' num2str(length(sessions))]);
        count = count + 1;
    end
end


%%
clear all
load /mnt/storage/rrr_magnum/M2/review_brain_states.mat

p(:, :, [20, 41, 47, 59]) = []; % outlier

x = repmat(f(:), [1, size(p, 2), size(p, 3)]);
c = repmat({'REST1', 'RUN', 'REST2'}, [size(p, 1), 1, size(p, 3)]);

g = gramm('x', x(:), 'y', p(:), 'color', c(:));
g.stat_summary('type', 'sem', 'geom', 'area');
g.geom_polygon('x', {[.05, .5, .5, .05]}, 'y', {[10^-3, 10^-3, 10^-1, 10^-1]});
% g.geom_polygon('x', {[1, 10, 10, 1]}, 'y', {[3*10^-6, 3*10^-6, 3*10^-4, 3*10^-4]});
g.axe_property('XScale', 'log', 'YScale', 'log', 'YLim', [10^-3, 10^-1], 'XLim', [f(1), f(end)]);
g.set_names('x', 'frequency (Hz)', 'y', 'PSD (Hz^2/Hz)', 'color', 'epoch');
g.set_text_options('base_size', 18);
figure
g.draw;

auc_so = trapz(f(f >= .05 & f <= .5), p(f >= .05 & f <= .5, :, :));
auc_fo = trapz(f(f >= 1), p(f >= 1, :, :));
auc_so = squeeze(auc_so)';
auc_fo = squeeze(auc_fo)';

p_so = zeros(3);
p_fo = zeros(3);
for ii = 1:3
    for jj = 1:3
        p_so(ii, jj) = signrank(auc_so(:, ii), auc_so(:, jj));
        p_fo(ii, jj) = signrank(auc_fo(:, ii), auc_fo(:, jj));
    end
end
idx = triu(true(length(p_so)), 1);
p_so = p_so(idx);
p_fo = p_fo(idx);

p_so = p_so * length(p_so); % bonferroni correction
p_so(p_so > 1) = 1;
p_fo = p_fo * length(p_fo); % bonferroni correction
p_fo(p_fo > 1) = 1;

figure
subplot(1, 2, 1)
boxplot(auc_so)
ylim([0 12*10^-3])
subplot(1, 2, 2)
boxplot(auc_fo)
ylim([0 5*10^-2])

%%
clear all
load /mnt/storage/rrr_magnum/M2/review_brain_states.mat

% p(:, :, [20, 41, 47, 59]) = []; % outlier

x = repmat(f(:), [1, size(coh, 2)]);

g = gramm('x', x(:), 'y', coh(:));
g.stat_summary('type', 'sem', 'geom', 'area');
% g.geom_polygon('x', {[.05, .5, .5, .05]}, 'y', {[.2, .2, .7, .7]});
% g.geom_polygon('x', {[1, 10, 10, 1]}, 'y', {[3*10^-6, 3*10^-6, 3*10^-4, 3*10^-4]});
g.axe_property('XScale', 'log', 'YLim', [.25 .65], 'XLim', [f(1), f(end)]);
g.set_names('x', 'frequency (Hz)', 'y', 'coherence');
g.set_text_options('base_size', 18);
figure
g.draw;

%%
clear all
load /mnt/storage/rrr_magnum/M2/review_brain_states.mat

p_lin(:, :, [20, 41, 47, 59]) = []; % outlier

x = repmat(f_lin(:), [1, size(p_lin, 2), size(p_lin, 3)]);
c = repmat({'REST1', 'RUN', 'REST2'}, [size(p_lin, 1), 1, size(p_lin, 3)]);

g = gramm('x', x(:), 'y', p_lin(:), 'color', c(:));
g.stat_summary('type', 'sem', 'geom', 'area');
g.geom_polygon('x', {[.05, .5, .5, .05]}, 'y', {[10^-3, 10^-3, 10^-1, 10^-1]});
% g.geom_polygon('x', {[1, 10, 10, 1]}, 'y', {[3*10^-6, 3*10^-6, 3*10^-4, 3*10^-4]});
% g.axe_property('XScale', 'log', 'YScale', 'log', 'YLim', [10^-3, 10^-1], 'XLim', [f(1), f(end)]);
g.set_names('x', 'frequency (Hz)', 'y', 'PSD (Hz^2/Hz)', 'color', 'epoch');
g.set_text_options('base_size', 18);
figure
g.draw;

auc = trapz(f_lin(f_lin <= 1), p_lin(f_lin <= 1, :, :));
auc = squeeze(auc)';

p_so = zeros(3);
for ii = 1:3
    for jj = 1:3
        p_so(ii, jj) = signrank(auc(:, ii), auc(:, jj));
    end
end
idx = triu(true(length(p_so)), 1);
p_so = p_so(idx);

p_so = p_so * length(p_so); % bonferroni correction
p_so(p_so > 1) = 1;

figure
boxplot(auc)
% ylim([0 12*10^-3])

%%
clear all
load /mnt/storage/rrr_magnum/M2/review_brain_states.mat

p_neur(:, :, [20, 41, 47, 59]) = []; % outlier

x = repmat(f(:), [1, size(p_neur, 2), size(p_neur, 3)]);
c = repmat({'REST1', 'RUN', 'REST2'}, [size(p_neur, 1), 1, size(p_neur, 3)]);

g = gramm('x', x(:), 'y', p_neur(:), 'color', c(:));
g.stat_summary('type', 'sem', 'geom', 'area');
g.geom_polygon('x', {[.05, .5, .5, .05]}, 'y', {[10^-3, 10^-3, 10^-1, 10^-1]});
% g.geom_polygon('x', {[1, 10, 10, 1]}, 'y', {[3*10^-6, 3*10^-6, 3*10^-4, 3*10^-4]});
% g.axe_property('XScale', 'log', 'YScale', 'log', 'XLim', [f(1), f(end)]);
% g.axe_property('XScale', 'log', 'YScale', 'log', 'YLim', [10^-3, 10^-1], 'XLim', [f(1), f(end)]);
g.set_names('x', 'frequency (Hz)', 'y', 'PSD (Hz^2/Hz)', 'color', 'epoch');
g.set_text_options('base_size', 18);
figure
g.draw;


%%
clear all

load 2/plane1/deconv.mat
d1 = deconv;
load 3/plane1/deconv.mat
d2 = deconv;

f = logspace(-2, 1, 2^8);

% fs = 1 / median(diff(tcs.tt));
fs = 19.01;

d1 = mean(d1, 2);
d2 = mean(d2, 2);

wdw = round(fs*200);
nover = round(wdw * .75);
p1 = pwelch(d1, hann(wdw), nover, f, fs, 'psd');
p2 = pwelch(d2, hann(wdw), nover, f, fs, 'psd');

figure
loglog(f, p1)
hold on
loglog(f, p2)


%%
clear all
load /mnt/storage/rrr_magnum/M2/review_brain_states.mat

rrr = zeros(length(R), 3);

for c = 1:length(R)
    idx = triu(true(size(R{c}, 1)), 1);
    for ii = 1:size(R{c}, 3)
        r = R{c}(:, :, ii);
        r = r(idx);
        
        rrr(c, ii) = mean(r);
    end
end

pfr = pfr';

p_pfr = zeros(3);
p_rrr = zeros(3);
for ii = 1:3
    for jj = 1:3
        p_pfr(ii, jj) = signrank(pfr(:, ii), pfr(:, jj));
        p_rrr(ii, jj) = signrank(rrr(:, ii), rrr(:, jj));
    end
end
idx = triu(true(length(p_pfr)), 1);
p_pfr = p_pfr(idx);
p_rrr = p_rrr(idx);

p_pfr = p_pfr * length(p_pfr); % bonferroni correction
p_pfr(p_pfr > 1) = 1;
p_rrr = p_rrr * length(p_rrr); % bonferroni correction
p_rrr(p_rrr > 1) = 1;

p_pfr
p_rrr

figure
subplot(1, 2, 1)
boxplot(pfr)
ylim([.2 1.2])
subplot(1, 2, 2)
boxplot(rrr)
ylim([0 .02])


%%
clear all

run = LFP('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_2.abf');
run.perform_analysis;

rest1 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_1.abf');
rest1.import_analysis(run.analysis)

rest2 = ensemble('/mnt/storage/HaoRan/RRR_motor/M2/RSC037/2017_08_18/2017_08_18_3.abf');
rest2.import_analysis(run.analysis)

%%
d1 = rest1.twop.deconv(:, run.analysis.order);
d1 = ca_filt(d1);
d1 = d1 > 0;

d2 = run.twop.deconv(:, run.analysis.order);
d2 = ca_filt(d2);
d2 = d2 > 0;

d3 = rest2.twop.deconv(:, run.analysis.order);
d3 = ca_filt(d3);
d3 = d3 > 0;

mua = cat(2, mean(d1, 2), mean(d2, 2), mean(d3, 2)) .* run.twop.fs;
mua = movmean(mua, run.twop.fs * 1);

tcs.tt = rest1.behavior.ts;
beh1 = convert_behavior(rest1.behavior, tcs, rest1.twop.deconv);

idx = discretize(rest1.behavior.ts, rest1.twop.ts);
vel1 = accumarray(idx(~isnan(idx))', rest1.behavior.speed(~isnan(idx)), [length(rest1.twop.ts), 1], @mean);
idx = discretize(run.behavior.ts, run.twop.ts);
vel2 = accumarray(idx(~isnan(idx))', run.behavior.speed(~isnan(idx)), [length(run.twop.ts), 1], @mean);
idx = discretize(rest2.behavior.ts, rest2.twop.ts);
vel3 = accumarray(idx(~isnan(idx))', rest2.behavior.speed(~isnan(idx)), [length(rest2.twop.ts), 1], @mean);
vel = cat(2, vel1, vel2, vel3);

colors = [255,   0,  7
           94, 169, 182
          105, 255, 75];
colors = colors ./ 255;

% d1 = repmat(d1, [1, 1, 3]) .* permute(colors(:, 1), [2 3 1]);
% d2 = repmat(d2, [1, 1, 3]) .* permute(colors(:, 2), [2 3 1]);
% d3 = repmat(d3, [1, 1, 3]) .* permute(colors(:, 3), [2 3 1]);

d = cat(4, d1, d2, d3);
d = repmat(d, [1, 1, 3, 1]);
d = d .* permute(colors, [3 4 1 2]);
d = permute(d, [2 1 3 4]);
for ii = 1:size(d, 4)
    mask = ~any(d(:, :, :, ii), 3);
    for jj = 1:size(d, 3)
        temp = d(:, :, jj, ii);
        temp(mask) = 1;
        d(:, :, jj, ii) = temp;
    end
end

lims = [6e3; 5e3; 1e3];
lims = cat(2, lims, lims + round(run.twop.fs * 60));

figure
for ii = 1:3
    ax(1) = subplot(5, 3, [1 4 7] + ii - 1);
    image(d(:, :, :, ii));
    ax(2) = subplot(5, 3, 10 + ii - 1);
    plot(mua(:, ii));
    ylim([0 1]);
    ax(3) = subplot(5, 3, 13 + ii - 1);
    plot(vel(:, ii));
    ylim([0 50]);
    
    linkaxes(ax, 'x');
    xlim(lims(ii, :));
end




