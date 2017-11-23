function [ev,rev] = explained_variance(rest1,run,rest2,sig,plotFlag)
% input deconv series (taking care to remove run epochs) and output ev, rev
rest1=fast_smooth(rest1,sig);
run=fast_smooth(run,sig);
rest2=fast_smooth(rest2,sig);

cor1=corr(rest1);
cor2=corr(run);
cor3=corr(rest2);

corr1=[];corr2=[];corr3=[];
for i=2:1:size(cor1,1)-1
    corr1=[corr1 cor1(i,1:i-1)];
    corr2=[corr2 cor2(i,1:i-1)];
    corr3=[corr3 cor3(i,1:i-1)];
end

[ev,p1]=partialcorr(corr1',corr2',corr3');
[rev,p2]=partialcorr(corr2',corr3',corr1');

if nargin>4 && plotFlag
    figure;
    subplot(1,2,1);
    plot(corr1,corr2,'b.');
    hold on
    plot(0:1,ev.*[0:1],'r--','linewidth',2);
    axis square
    xlabel('Run Correlations')
    ylabel('Rest 1 (pre-exposure) Correlations')
    title(['Corr Coef = ' num2str(ev) ' REV = ' num2str(ev^2) ' p = ' num2str(p1)]);
    xlim([-0.2 0.8]);
    ylim([-0.2 0.8]);


    subplot(1,2,2);
    plot(corr3,corr2,'b.');
    hold on
    plot(0:1,rev.*[0:1],'r--','linewidth',2);
    axis square
    xlabel('Run Correlations')
    ylabel('Rest 2 (post-exposure) Correlations')
    title(['Corr Coef = ' num2str(rev) ' EV = ' num2str(rev^2) ' p = ' num2str(p2)]);
    xlim([-0.2 0.8]);
    ylim([-0.2 0.8]);
end

ev=ev^2;
rev=rev^2;