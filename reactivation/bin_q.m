function q = bin_q(deconv,win)
% bin Q-matrix with window size 'win'

q=zeros(size(deconv,1)/win,size(deconv,2));

count=0;
while(count<size(deconv,1)/win)
    q(count+1,:)=sum(deconv(count*win+1:(count+1)*win,:));
    count=count+1;
end

try
    q=[q;sum(deconv(count*win+1:end,:))];
catch
end

q=q>0;