function [r, lags] = bxcorr2(A, B, maxlag, scaleopt, numThreads)
% by default, stupid matlab xcorr2 doesn't allow definition of maxlag
% this could result in reduced performance for certain applications
% this little script works exactly like xcorr, but for the 2D case
% maxlag can be a scalar, or a two-elements vector
% scaleopt: 'none', 'biased', 'unbiased', 'normalized'/'coeff'
% numThreads: number of threads to use (works well with hyperthreading)
%
% matxcorr2.cpp needs to be compiled before usage.

if nargin < 3
    maxlag = floor( ( sum( [ size(A); size(B) ] ) - 1 ) ./ 2 );
end
if isempty(maxlag)
    maxlag = floor( ( sum( [ size(A); size(B) ] ) - 1 ) ./ 2 );
end
if nargin < 4
    scaleopt = 'none';
end
if nargin < 5
    numThreads = maxNumCompThreads*2;
    if isunix
        [~,numThreads] = system('grep -c ^processor /proc/cpuinfo');
        numThreads = str2double(numThreads);
%         !uname -or
%         !echo "using" $(grep -c ^processor /proc/cpuinfo) "logical processing cores"
    elseif ispc
        disp('You should switch to Linux you Windows peasant!')
    end
end
if size(maxlag) < 2
    maxlag = ones(1,2) .* maxlag;
end

lags.x = -maxlag(1):maxlag(1);
lags.y = -maxlag(2):maxlag(2);
if mod(size(A,1) + size(B,1), 2)
    lags.x(lags.x==0) = [];
end
if mod(size(A,2) + size(B,2), 2)
    lags.y(lags.y==0) = [];
end

r = matxcorr2(A, B, int16(maxlag), numThreads);

switch lower(scaleopt)
    case 'none'
        return;
    case 'biased'
        r = r ./ matxcorr2(ones(size(A)), ones(size(B)), int16([0 0]), numThreads);
    case 'unbiased'
        r = r ./ matxcorr2(ones(size(A)), ones(size(B)), int16(maxlag), numThreads);
    case {'normalized','coeff'}
        rxx = matxcorr2(A, A, int16([0 0]), numThreads);
        ryy = matxcorr2(B, B, int16([0 0]), numThreads);
        r = r ./ sqrt(rxx * ryy);
    otherwise
        error(['scaleopt ''' scaleopt ''' is undefined']);
end