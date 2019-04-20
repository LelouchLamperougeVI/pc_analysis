function smoothed=fast_smooth(data,sig,dim)
% Faster 1d gaussian kernel smoothing
% Accounts for edge underestimation and NaNs
% Usage: 
%   data: matrix of size n x m where rows contain observations and columns
%      contain variables (can also be vector)
%   sig: kernel standard dev
%   dim: dimension of observations (default 1)
%
%   smoothed: smoothed data
%
% by HaoRan Chang

if sig==0
    smoothed=data;
    return
end

if nargin==2
    dim=1;
end

switch dim
    case 1
    case 2
        data=data';
    otherwise
        error('WTF dude?')
end

dataSize=size(data,1);
kernelSize=ceil(10*sig);
kernelSize=kernelSize-~mod(kernelSize,2); %avoid even kernel sizes
alpha=(kernelSize-1)/sig/2;
kernel=gausswin(kernelSize,alpha);
kernel=kernel./sum(kernel);

nan_idx = any(isnan(data),2);
taper=zeros(kernelSize,1);
data=[repmat(taper,1,size(data,2));data;repmat(taper,1,size(data,2))];
data=data(:);
data(isnan(data))=0;

smoothed=conv(data,kernel,'same');
smoothed=reshape(smoothed,dataSize+kernelSize*2,[]);
smoothed=smoothed./conv([zeros(kernelSize,1); ~nan_idx; zeros(kernelSize,1)],kernel,'same'); %account for edge data underestimation
smoothed([1:kernelSize end-kernelSize+1:end],:)=[];
smoothed(nan_idx,:) = nan;

switch dim
    case 1
    case 2
        smoothed=smoothed';
end