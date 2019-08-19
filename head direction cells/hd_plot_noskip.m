function hd_plot_noskip(hObj,evt)
    global n
    global h
    global trials_ts
    global deconv
    global den
    global unit_rot
    global mean_vect_rho
    global mean_vect_theta
    global pval
    global fn
    
    bin=100;
    
    switch evt.Key
        case 'rightarrow'
            n=n+1;
            if n>size(deconv,2)
                n=size(deconv,2);
            end
        case 'leftarrow'
            n=n-1;
            if n<1
                n=1;
            end
        otherwise
            return
    end
    figure(h);
    
    stack=zeros(bin,length(trials_ts)-1);
    
%     stack=zeros(max(diff(trials_ts))+1,length(trials_ts)-1);
    count=1;
    idx=zeros(1,length(trials_ts)-1);
%     for i=1:length(trials_ts)-1
    theta=0:360/bin:360;
    for i=1:length(trials_ts)-1
        signal=smooth(den(trials_ts(i):trials_ts(i+1),n),21);
        stack(:,i)=arrayfun(@(x) mean(signal(unit_rot(trials_ts(i):trials_ts(i+1))>=theta(x) & unit_rot(trials_ts(i):trials_ts(i+1))<theta(x+1))),1:bin);
    end
    stack=flip(stack);
    stack=rot90(stack,-1);
    stack(isnan(stack))=0;

    subplot(2,1,1);
    imagesc(stack);
    colorbar;
    set(gca,'xtick',0:size(stack,2)/6:size(stack,2));
    set(gca,'xticklabel',{'0','60','120','180','240','300','360'});
    xlabel('degrees'); ylabel('trial');
    
    subplot(2,1,2);
    polarplot(deg2rad(theta(1:end-1)),sum(stack)./size(stack,2)); hold on;
    polarplot([0 mean_vect_theta(1,n)],[0 mean_vect_rho(1,n)]); hold off;
    title(['p = ' mat2str(pval(n))]);
    
    supertitle({fn, ['n = ' mat2str(n)]});
    colormap hot;