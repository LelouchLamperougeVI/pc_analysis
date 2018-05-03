function H=loc_entropy(deconv,win,precision)
% Calculate local entropy within neighboring window
% Calculated at -log2(1/(precision+1)) bits precision
% Entropy normalized to 1 bit as max entropy

if nargin<3
    precision=5;
end
win=floor(win/2);

deconv=(deconv-min(deconv))./max(deconv-min(deconv));
deconv=round(deconv.*precision);

% H=zeros(size(deconv,1),1);
p=zeros(size(deconv,1),precision+1);
for i=0:precision
    p(:,i+1)=sum(conv2(deconv==i,ones(win*2+1,1),'same'),2);
end
p=p./sum(p,2);
p([1:win end-win+1:end],:)=0;
H=-p.*log2(p);
H(isnan(H))=0;
H=sum(H,2);
% H=H./(-log2(1/(precision+1)));