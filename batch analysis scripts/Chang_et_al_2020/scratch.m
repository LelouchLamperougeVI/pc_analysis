span = 1;
wdw = round(lfp.lfp.fs .* span);

swr_specs = [];
for sw_id = 1:length(lfp.lfp.swr.swr_on)
    idx = find(lfp.lfp.ts == lfp.lfp.swr.swr_on(sw_id), 1);
    try
        signal = lfp.lfp.lfp(idx - wdw : idx + wdw);
    catch
        continue
    end
    [s, f, coi] = cwt(signal, 'amor', lfp.lfp.fs);
    swr_specs = cat(3, swr_specs, s);
    disp([num2str(sw_id) '/' num2str(length(lfp.lfp.swr.swr_on))]);
end

null = [];
for ii = 1:length(lfp.lfp.swr.swr_on)
    on = randperm(length(lfp.lfp.lfp) - length(signal), 1);
    s = cwt(lfp.lfp.lfp(on : on + length(signal) - 1), 'amor', lfp.lfp.fs);
    null = cat(3, null, s);
    disp([num2str(ii) '/' num2str(length(lfp.lfp.swr.swr_on))]);
end


figure
% spec = abs( (mean(log(swr_specs), 3) - mean(log(null), 3)) ./ std(log(null), [], 3) );
spec = abs( mean(swr_specs, 3) );
% spec = abs( (mean((swr_specs), 3) - mean((null), 3)) ./ std((null), [], 3) );
t = linspace(-span, span, size(spec, 2));
% imagesc('xdata', t, 'ydata', f(f<=300 & f>=100), 'cdata', spec(f<=300 & f>=100, :));
imagesc('xdata', t, 'ydata', f(f<=300), 'cdata', spec(f<=300, :));
% imagesc('xdata', t, 'ydata', f(f<=250), 'cdata', abs(mean(log(swr_specs(f<=250, :, :)), 3)));
xline(0);
colormap jet


%% SWR FIR filter
n = 400; % filter order
cF = [150 250]; % cutoff frequencies
wdw = hamming(n+1); % creat hamming window
flt = fir1(n, cF/(lfp.lfp.fs/2), 'bandpass', wdw, 'scale');
disp(['Filter is linear-phase: ' num2str(islinphase(flt))])

freqz(flt, 1, 512);

test = filtfilt(flt, 1, lfp.lfp.lfp);

figure
plot(lfp.lfp.ts, lfp.lfp.lfp);
hold on
plot(lfp.lfp.ts, test);

%% SWR envelope

y = lfp.lfp.swr.swr;
ts = lfp.lfp.ts;
wdw = round(lfp.lfp.fs * .008); % 8 ms RMS window

env = envelope(y, wdw, 'rms');

swr_thres = 3;
swr_idx = env > (mean(env, 'omitnan') + swr_thres*std(env, 'omitnan'));
onset = env > (mean(env, 'omitnan') + .75*swr_thres*std(env, 'omitnan'));

figure;
plot(ts, y);
hold on
plot(ts, env);
plot(ts(get_head(swr_idx)), ones(sum(get_head(swr_idx)), 1) .* max(env), 'k*')



