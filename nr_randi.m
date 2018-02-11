function X=nr_randi(imax,m,n)
  % Generate uniformly distributed pseudorandom integers from 1 to imax without replacement
  if imax<m*n
    error('sample exceeded data range');
  end

  X=zeros(m,n);
  sample=1:imax;
  X=datasample(sample,m*n,'Replace',false);
  X=reshape(X,m,n);
