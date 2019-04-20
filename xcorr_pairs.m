function [r, lags] = xcorr_pairs(A, maxlag)

r = zeros(2*maxlag+1, size(A,2) * (size(A,2) - 1) / 2);

count = 1;
for i = 1:size(A, 2)
    for j = (1+i):size(A,2)
        r(:, count) = xcorr(A(:,i), A(:,j), maxlag, 'coeff');
        count = count+1;
    end
end

lags = -maxlag:maxlag;