%% belt
bins = 50;
belt_length = 144.5;
% cue_centres = mod([0 15.5 30 60 115] + (belt_length - 30), belt_length);
cue_centres = mod([0 15.5 60 115] + (belt_length - 30), belt_length);
% cue_widths = [2.5 5 1 5.5 3];
cue_widths = [2.5 5 5.5 3];

belt_idx = false(1, bins);

for i = 1:length(cue_centres)
    idx = [cue_centres(i)-cue_widths(i)/2 cue_centres(i)+cue_widths(i)/2];
    idx = knnsearch(linspace(0, belt_length, bins)', idx');
    belt_idx(idx(1):idx(2)) = true;
end