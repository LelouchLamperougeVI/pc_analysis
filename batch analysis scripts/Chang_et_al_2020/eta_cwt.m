function [s, t, f, h] = eta_cwt(varargin)
% Event-triggered average Continuous Wavelet Tranform
% Compute CWT -> zscore wavelet magnitude -> obtain average zscored TF
% Throws away edges for acurate estimation
%
% Inputs
%   y:      time-series signal
%   fs:     signal sampling rate (in Hz)
%   wdw:    window for observation of events [-wdw, wdw] (in sec)
%   evts:   event onset time stamps (in sec)
%
% Name, Value pairs
%   'nans': logical vector of rejection epophs (same length as y)
%   'plot': plot resulting spectrogram (default false)
%   'nvc':  number of voices per octave (default 10)
%
% Outputs
%   s:      spectrogram of z-scored wavelet magnitude
%   t:      time-stamps of spectrogram
%   f:      frequencies of spectrogram
%   h:      figure handle (if 
%
% by HaoRan Chang, PhD candidate
% Polaris Research Group, Canadian Centre for Behavioural Neuroscience
% University of Lethbridge, Alberta, Canada

[y, fs, wdw, evts, nans, nvc, plt] = parse_inputs(varargin);

[s, f, coi] = cwt(y, 'amor', fs, 'voicesperoctave', nvc);
t = linspace(0, (length(y) - 1) / fs, size(s, 2))';

% remove edge artifacts
nan_coi = inf(length(coi), 1);
heads = find(get_head(~nans));
tails = get_head(~nans(end:-1:1));
tails = find(tails(end:-1:1));

for ii = 1:length(heads)
    span = heads(ii) : heads(ii) + floor((tails(ii) - heads(ii) ) / 2);
    nan_coi(span) = coi(1:length(span));
    span = heads(ii) + floor((tails(ii) - heads(ii) ) / 2) + 1 : tails(ii);
    nan_coi(span) = coi(end - length(span) + 1:end);
end

coi = double(nan_coi' < f);
coi(~coi) = nan;
s = abs(s) .* coi;

zscored = ( s - mean(s, 2, 'omitnan') ) ./ std(s, [], 2, 'omitnan');
wdw = round(wdw / median(diff(t)));

mean_zscore_spec = zeros(size(s,1), wdw*2 + 1, length(evts));
idx = knnsearch(t, evts);
for ii = 1:length(evts)
    mean_zscore_spec(:,:,ii) = zscored(:, idx(ii) - wdw : idx(ii) + wdw);
end
mean_zscore_spec = mean(mean_zscore_spec, 3, 'omitnan');

t = linspace(-median(diff(t)) * wdw, median(diff(t)) * wdw, wdw*2 + 1);
s = mean_zscore_spec;

if plt
    h = figure;
    [X, Y] = meshgrid(t, f);
    ph = pcolor(X, Y, s);
    ph.EdgeColor = 'none';
    ph.FaceColor = 'interp';
    colormap jet
    c = colorbar;
    c.Label.String = 'z-score Wavelet Magnitude';
    xline(0, 'linewidth', 2, 'color', [1 1 1])
    xlabel('time from onset (sec)');
    ylabel('frequency (Hz)');
else
    h = [];
end


function [y, fs, wdw, evts, nans, nvc, plt] = parse_inputs(inputs)
if length(inputs) < 4
    error('Minimum of 4 inputs required.');
end
y = inputs{1};
fs = inputs{2};
wdw = inputs{3};
evts = inputs{4};
nans = false(length(y), 1);
nvc = 10;
plt = false;

count = 5;
while(count < length(inputs))
    switch lower(inputs{count})
        case 'nans'
            nans = inputs{count+1};
        case 'nvc'
            nvc = inputs{count+1};
        case 'plot'
            plt = inputs{count+1};
        otherwise
            error([inputs{count} ' is not a valid parameter.']);
    end
    count = count + 2;
end