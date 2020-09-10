span = .3;
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
% spec = abs( mean(swr_specs, 3) );
% spec = abs( (mean((swr_specs), 3) - mean((null), 3)) ./ std((null), [], 3) );
spec = abs( (mean((swr_specs), 3) - mean(null, [2 3])) ./ std((null), [], [2 3]) );
t = linspace(-span, span, size(spec, 2));
% imagesc('xdata', t, 'ydata', f(f<=300 & f>=100), 'cdata', spec(f<=300 & f>=100, :));
imagesc('xdata', t, 'ydata', f(f<=300), 'cdata', spec(f<=300, :));
% imagesc('xdata', t, 'ydata', f(f<=250), 'cdata', abs(mean(log(swr_specs(f<=250, :, :)), 3)));
xline(0);
colormap jet



%%
[s, t, f] = eta_cwt(lfp.lfp.lfp, lfp.lfp.fs, .2, lfp.lfp.swr.swr_on, 'nans', isnan(lfp.lfp.swr.swr_env), 'plot', true, 'nvc', 20);
[X, Y] = meshgrid(t, f);

figure
h = pcolor(X, Y, s);
h.EdgeColor = 'none';
h.FaceColor = 'interp';
ylim([0 300]);
caxis([0 4]);
colormap jet
c = colorbar;
c.Label.String = 'z-score Wavelet Magnitude';
xline(0, 'linewidth', 2, 'color', [1 1 1])



