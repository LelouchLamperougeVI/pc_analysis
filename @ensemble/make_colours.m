function make_colours(obj)
% Any colour you like - can't really think of a description for this
% function so there's a Pink Floyd reference for you instead :)
% It's a good song!

if isempty(obj.ensembles.clust)
    error('you need to cluster the data first...');
end

obj.ensembles.colours = distinguishable_colors(length(obj.ensembles.clust), {'w','k'});