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
    [~,~,x]=par_sort(x);
    [~,~,y]=par_sort(y);
end

rho=corr(x,y);
