clear all
list = {'RSC036', 'RSC037', 'RSC038'};
bias = [];
for f = 1:length(list)
    root = dir(list{f});
    root(1:2) = [];
    
    for s = 1:length(root)
        lfp = LFP(fullfile(root(s).folder, root(s).name, [root(s).name '_2.abf']));
        try
            review_round2;
        catch
        end
    end
end