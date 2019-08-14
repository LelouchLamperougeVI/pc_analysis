function [r, lags] = bxcorr2(A, B, maxlag, scaleopt, numThreads)
% by default, stupid matlab xcorr2 doesn't allow definition of maxlag
% this could result in reduced performance for certain applications
% this little script works exactly like xcorr, but for the 2D case
% maxlag can be a scalar, or a two-elements vector
% scaleopt: 'none', 'biased', 'unbiased', 'normalized'/'coeff'
% numThreads: number of threads to use (works well with hyperthreading)
%
% matxcorr2.cpp needs to be compiled before usage.

dims = min( [ size(A); size(B) ]);
if nargin < 3
    maxlag = dims - 1;
end
if isempty(maxlag)
    maxlag = dims - 1;
end
if nargin < 4
    scaleopt = 'none';
end
if nargin < 5
    numThreads = maxNumCompThreads*2;
end
if size(maxlag) < 2
    maxlag = ones(1,2) .* maxlag;
end

lags.x = -maxlag(1) : maxlag(1);
lags.y = -maxlag(2) : maxlag(2);

if mod(size(A,1) + size(B,1), 2)
    lags.x(lags.x==0)=[];
    dims(1) = dims(1) + 1;
end
if mod(size(A,2) + size(B,2), 2)
    lags.y(lags.y==0)=[];
    dims(2) = dims(2) + 1;
end

r = matxcorr2(A, B, int16(maxlag), numThreads);

switch lower(scaleopt)
    case 'none'
        return;
    case 'biased'
        r = r ./ (numel(A) + numel(B));
    case 'unbiased'
        r = r ./ ( (dims(1) - abs(lags.x)') * (dims(2) - abs(lags.y)) );
    case {'normalized','coeff'}
        
    otherwise
        error(['scaleopt ''' scaleopt ''' is undefined']);
end