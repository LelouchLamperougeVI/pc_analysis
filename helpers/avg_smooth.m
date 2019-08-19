function smoothed=avg_smooth(data,k,dim)
% smooth each column using averaging kernel of size k
% NaNs are omitted

if k==0
    smoothed=data;
    return
end

if nargin==2
    dim=1;
end

if nargin<2
    error('Need 2 inputs')
end

if ~mod(k,2)
    warning('Unless circumstance dictates, it''s advised use an odd kernel')
end

switch dim
    case 1
    case 2
        data=data';
    otherwise
        error('Sorry. Only up to 2 dimensions are supported for now.')
end

idx=min(data)-1;
data=data-idx;
data(isnan(data))=0;

kernel=ones(k,size(data,2));
smoothed=mat_conv(data,kernel);

n=mat_conv(double(data>0),kernel);

smoothed=smoothed./n;
smoothed(smoothed==0)=nan;
smoothed=smoothed+idx;

smoothed([1:floor(k/2)-~mod(k,2) end-floor(k/2)+1:end],:)=[];

switch dim
    case 1
    case 2
        smoothed=smoothed';
end