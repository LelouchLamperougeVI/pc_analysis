function X = get_middle(X)

X(1,:)=0;
X(end,:)=0;

idx=diff(X);
idx=ceil((find(idx==-1)-find(idx==1))./2+find(idx==1));
X(1:end-1,:)=0;
X(idx)=1;
