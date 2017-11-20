cor1=corr(c1);
cor2=corr(c2);
cor3=corr(c3);

corr1=[];corr2=[];corr3=[];
for i=2:1:size(cor1,1)-1
    corr1=[corr1 cor1(i,1:i-1)];
    corr2=[corr2 cor2(i,1:i-1)];
    corr3=[corr3 cor3(i,1:i-1)];
end

[PreEx,p1]=partialcorr(corr1',corr2',corr3');
[PostEx,p2]=partialcorr(corr2',corr3',corr1');

subplot(1,2,1);
plot(corr1,corr2,'b.');
hold on
plot(0:1,PreEx.*[0:1],'r--','linewidth',2);
axis square
xlabel('Run Correlations')
ylabel('Rest 1 (pre-exposure) Correlations')
title(['Corr Coef = ' num2str(PreEx) ' REV = ' num2str(PreEx^2) ' p = ' num2str(p1)]);
xlim([-0.2 0.8]);
ylim([-0.2 0.8]);


subplot(1,2,2);
plot(corr3,corr2,'b.');
hold on
plot(0:1,PostEx.*[0:1],'r--','linewidth',2);
axis square
xlabel('Run Correlations')
ylabel('Rest 2 (post-exposure) Correlations')
title(['Corr Coef = ' num2str(PostEx) ' EV = ' num2str(PostEx^2) ' p = ' num2str(p2)]);
xlim([-0.2 0.8]);
ylim([-0.2 0.8]);