function A = fill_gaps(A, distance)
% given a logical vector A, fill in gaps no longer than 'distance'

idx=find(A);

gaps=find(diff(idx)<=(distance+1) & diff(idx)>1);

for ii=1:length(gaps)
    A(idx(gaps(ii)):idx(gaps(ii)+1))=true;
end