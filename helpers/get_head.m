function X = get_head(X)
% Find the first element in a series of 1's

  idx=diff(X);
  idx=idx>0;
  X=[X(1,:);idx];
