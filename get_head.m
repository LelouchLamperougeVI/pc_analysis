function X = get_head(X)

  idx=diff(X);
  idx=idx>0;
  X=[X(1,:);idx];
