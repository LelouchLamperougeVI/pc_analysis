% function seq_align(temp,match)
sig=.2; %seconds
shuffles=1;

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

d=fast_smooth(deconv,19.1*sig);
d=zscore(d);
[coef,score]=pca(d);
cor=zeros(size(d,2),6);
for i=1:6
    cor(:,i)=arrayfun(@(x) xcorr(d(:,x),score(:,i),0), 1:size(d,2));
end

%%
ax1=subplot(2,1,1);
imagesc(d');
ylabel('neuron no.')
ax2=subplot(2,1,2);
plot(score(:,1:6));
xlabel('frame')
ylabel('PC score')
linkaxes([ax1 ax2],'x');