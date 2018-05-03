function H=bin_loc_entropy(deconv,win)

win=floor(win/2);

num_microstates=min([2^(win*2+1) size(deconv,2)]);

H=zeros(size(deconv,1),1);
for i=win+1:size(deconv,1)-win
    temp=deconv(i-win:i+win,:);
    temp=bi2de(temp);
    
    p=unique(temp);
    
    p=repmat(temp',1,length(p))==p;
    p=sum(p)./size(deconv,2);
    
    H(i)=-sum(p.*log2(p));
    
%     if H(i)>log2(num_microstates)
%         error('wtf? entropy exceeds maximum entropy of system')
%     end
end

% H=H./log2(num_microstates);