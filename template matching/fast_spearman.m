function rho=fast_spearman(x,y,gpuFlag)
% A little bit faster than corr(x,y,'type','spearman')
% Expect significant performance increase with GPU

if nargin<3
    gpuFlag=false;
end

if gpuFlag
    A=(A-min(A))./range(A).*(2^16-2);
    A=uint16(A);
    [~,~,ranks]=oddeven(A);
    ranks=double(gather(ranks));
else
    [sorted,idx]=par_sort(x);
    ranks=get_ranks(sorted);
    idx=idx+(0:size(x,2)-1).*size(x,1);
    x(idx)=ranks;
    [sorted,idx]=par_sort(y);
    ranks=get_ranks(sorted);
    idx=idx+(0:size(y,2)-1).*size(y,1);
    y(idx)=ranks;
end

rho=corr(x,y);
