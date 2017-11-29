function smoothed=fast_smooth(data,sig,dim)
% Faster 1d gaussian kernel smoothing
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
alpha=(kernelSize-1)/sig/2;
kernel=gausswin(kernelSize,alpha);
kernel=kernel./sum(kernel);

taper=zeros(kernelSize*2,1);
data=[data;repmat(taper,1,size(data,2))];
data=reshape(data,size(data,1)*size(data,2),1);

smoothed=conv(kernel,data);
smoothed=smoothed(floor(kernelSize/2)+1:end-ceil(kernelSize/2)+1);
smoothed=reshape(smoothed,dataSize+kernelSize*2,[]);
smoothed(end-kernelSize*2+1:end,:)=[];

switch dim
    case 1
    case 2
        smoothed=smoothed';
end