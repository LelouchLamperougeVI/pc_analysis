function [scores,pval]=mult_needleWunsch(sequences)
% multiple pairwise comparisons of sequences
% generates score matrix of pairwise comparisons
%
% Inputs
%   sequences:  a N x M matrix with N individual sequences of length M
%
% Output
%   scores: matrix of Needleman-Wunsch scores for sequence pairs
%   pval:   significance of score calculated from shuffled dist

scores=needleWunsch(sequences(1,:),sequences(1,:));
scores=diag(scores.*ones(1,size(sequences,1)));

for i=1:size(sequences,1)
    for j=i+1:size(sequences,1)
        scores(i,j)=needleWunsch(sequences(i,:),sequences(j,:));
    end
end

scores=scores+scores'-diag(diag(scores));

shuffles=1000;
shuffled_scores=zeros(1,shuffles);
for i=1:shuffles
    shuffled_scores(i)=needleWunsch(sequences(1,randperm(size(sequences,2))),sequences(1,randperm(size(sequences,2))));
end

pval=zeros(size(sequences,1));
for i=1:size(pval,1)
    for j=1:size(pval,2)
        pval(i,j)=1-sum(scores(i,j)>shuffled_scores)/shuffles;
    end
end