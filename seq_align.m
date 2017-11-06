function seq_align(temp,match)
sig=5;
shuffles=10;

coef=zeros(size(deconv,2),size(deconv,2),shuffles+1);
coef(:,:,1)=pca(Smooth(deconv,[19.1*sig 0]));

for shuffles=1:shuffles
    shifted_d=deconv;
    for i=1:size(deconv,2)
        shift=ceil(size(deconv,1)*rand(1));
        shifted_d(:,i)=circshift(shifted_d(:,i),shift);
    end
    coef(:,:,shuffles+1)=pca(Smooth(shifted_d,[19.1*sig 0]),'algorithm','eig');
end