function [sequence,out_of_bound]=get_sequence(cells,interp_ts,t,ts,win)

q=interp_ts(:,cells);
idx1=find(t-win==ts);
if isempty(idx1)
    idx1=1;
end
idx2=find(t+win==ts);
if isempty(idx2)
    idx2=size(q,1);
end
q=q(idx1:idx2,:);

[~,sequence]=min(q);

out_of_bound=sequence==1 | sequence==idx2-idx1+1;

sequence=sequence+idx1-1;

% m=median(q,2);

% sequence=zeros(1,size(q,2));
% for i=1:size(q,2)
%     C=xcorr(m,q(:,i));
%     [~,sequence(i)]=max(C);
% end