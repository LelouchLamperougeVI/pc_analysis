function [ev,rev] = explained_variance(deconv,behavior,sig,plotFlag)
% input deconv series (taking care to remove deconv{2} epochs) and output ev, rev

idx={abs(behavior{1}.speed_raw)>noRun(behavior{2}.speed_raw),abs(behavior{3}.speed_raw)>noRun(behavior{2}.speed_raw)};
% sum(idx{1})
% sum(idx{2})

deconv{1}=fast_smooth(deconv{1},sig);
deconv{2}=fast_smooth(deconv{2},sig);
deconv{3}=fast_smooth(deconv{3},sig);

deconv{1}(idx{1},:)=[];
deconv{3}(idx{2},:)=[];

cor1=corr(deconv{1});
cor2=corr(deconv{2});
cor3=corr(deconv{3});

corr1=[];corr2=[];corr3=[];
for i=2:1:size(cor1,1)-1
    corr1=[corr1 cor1(i,1:i-1)];
    corr2=[corr2 cor2(i,1:i-1)];
    corr3=[corr3 cor3(i,1:i-1)];
end

[rev,p1]=partialcorr(corr1',corr2',corr3');
[ev,p2]=partialcorr(corr2',corr3',corr1');

if nargin>3 && plotFlag
    figure;
    subplot(1,2,1);
    plot(corr1,corr2,'b.');
    hold on
    plot(0:1,rev.*[0:1],'r--','linewidth',2);
    axis square
    xlabel('deconv{2} Correlations')
    ylabel('Rest 1 (pre-exposure) Correlations')
    title(['Corr Coef = ' num2str(rev) ' REV = ' num2str(rev^2) ' p = ' num2str(p1)]);
    xlim([-0.2 0.8]);
    ylim([-0.2 0.8]);


    subplot(1,2,2);
    plot(corr3,corr2,'b.');
    hold on
    plot(0:1,ev.*[0:1],'r--','linewidth',2);
    axis square
    xlabel('deconv{2} Correlations')
    ylabel('Rest 2 (post-exposure) Correlations')
    title(['Corr Coef = ' num2str(ev) ' EV = ' num2str(ev^2) ' p = ' num2str(p2)]);
    xlim([-0.2 0.8]);
    ylim([-0.2 0.8]);
end

ev=ev^2;
rev=rev^2;