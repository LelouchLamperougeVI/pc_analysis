function chan = detect_channels(obj)
% Automatically detect channels

d = obj.abf.raw;
si = obj.abf.si;
Fs = 1 / (si / 1e6);
s = downsample(d, floor(Fs/1e3));
Fs = size(s, 1) / size(d, 1) * Fs;
ndft = 2^12;
T = abs(fft(s, ndft));
f = linspace(0, Fs, size(T, 1));

twopFrange = [3 40]; % range of frequencies typical for 2p

bmratio = (sum(d>1) ./ sum(d<1));
steps = get_head(d>1);
freq = arrayfun(@(x) mean(diff(find(steps(:,x)))), 1:length(bmratio));
freq = 1 ./ ((si / 1e6) .* freq);
step_count = sum(steps);
[~, guess] = sort(bmratio, 'descend');

guess(range(s(:, guess)) < 1) = []; % silent channels;

if sum(bmratio > 10) > 1
    error('More than one camera/2photon channels exist. Automatic detection failed.');
elseif sum(bmratio > 10) < 1
    error('Failed to detect any 2p channel. I certainly hope you recorded the frames :P');
elseif mod(step_count(guess(1)), 10)
    error('Odd frame count for 2p.');
elseif freq(guess(1)) < twopFrange(1) || freq(guess(1)) > twopFrange(2)
    error('Unexpected 2p frequency.');
end
chan(1) = guess(1);
guess(1) = [];

if sum(bmratio(guess) > .1) ~= 2
    error('Failed to detect encoder channels. Make sure this is not resting data. Or perhaps you just got a really lazy animal :P');
end
a = guess(1);
b = guess(2);
guess(1:2) = [];
r = xcorr(s, round(Fs * .01), 'unbiased');
r = reshape(r, size(r,1), size(s,2), size(s,2));
r = (r - mean(r)) ./ range(r);
if sum( r(1 : floor(size(r,1)/2), a, b) ) > sum( r(end : -1 : ceil(size(r,1)/2) + 1, a, b) )
    temp = a;
    a = b;
    b = temp;
end
chan(2) = a;
chan(3) = b;

idx = abs(freq(guess) - step_count(guess) ./ (length(s) / Fs)) < 1e-3 & abs(median(s(:, guess))) < .1;
if sum(idx) == 0
    error('Failed to detect reward channel. Make sure this is not resting data. Or perhaps you just got a really lazy animal :P');
end
[~, rev] = min(bmratio(guess(idx)));
chan(4) = guess(rev);
guess(idx) = [];

f200 = f(f < 200 & f > 0);
md = fitlm(log(T(f < 200 & f > 0, guess)), log(f200), 'linear');
l = find(md.Coefficients.Estimate(2:end) < -0.1 & md.Coefficients.pValue(2:end) < 1e-6);

if isempty(l)
    warning('No LFP channel detected.');
    chan(6) = NaN;
else
    chan(6) = setdiff(guess(l), chan);
end

chan(5) = NaN;
chan(7) = NaN;
