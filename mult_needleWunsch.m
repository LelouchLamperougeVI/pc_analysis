function [scores,pval]=mult_needleWunsch(sequences)
% multiple pairwise comparisons of sequences
% generates score matrix of pairwise comparisons
%
% Inputs
%   sequences:  cell containing individual sequences
%
% Outputs
%   scores: matrix of Needleman-Wunsch scores for sequence pairs
%   pval:   significance of score calculated from shuffled dist
%
% by HaoRan Chang

% scores=needleWunsch(sequences{1},sequences{1});
% scores=diag(scores.*ones(1,length(sequences)));
scores=zeros(length(sequences));
for i=1:length(sequences)
    parfor j=i:length(sequences)
        scores(i,j)=needleWunsch(sequences{i},sequences{j});
    end
end

scores=scores+scores'-diag(diag(scores));

shuffles=1000;
shuffled_scores=zeros(1,shuffles);
parfor i=1:shuffles
    idx=[randi(length(sequences)) randi(length(sequences))];
    shuffled_scores(i)=needleWunsch(sequences{idx(1)}(randperm(length(sequences{idx(1)}))),sequences{idx(2)}(randperm(length(sequences{idx(2)}))));
end

pval=zeros(length(sequences));
for i=1:size(pval,1)
    for j=1:size(pval,2)
        pval(i,j)=1-sum(scores(i,j)>shuffled_scores)/shuffles;
    end
end