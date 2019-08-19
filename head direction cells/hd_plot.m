function hd_plot(hObj,evt)
    global n
    global h
    global frame_ts
    global trials_ts
    global trials
    global deconv
    global den
    global unit_rot
    global mean_vect_rho
    global mean_vect_theta
    global pval
    global fn
    
    bin=200;
    
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
%     for i=1:bin
%         try
%             stack(:,i)=smooth(deconv(trials_ts(i):trials_ts(i+1),n),21);
%         catch
%             stack(1:end-1,i)=smooth(deconv(trials_ts(i):trials_ts(i+1),n),21);
%         end
        if count<=length(trials)
            if frame_ts(trials_ts(i+1))>trials(count)
                count=count+1;
                idx(i)=1;
            end
        end
    end
    stack=flip(stack);
    stack=rot90(stack,-1);

    subplot(2,2,1);
    imagesc(stack(~idx,:));
    colorbar;
    title('No Reward');
    set(gca,'xtick',0:size(stack,2)/6:size(stack,2));
    set(gca,'xticklabel',{'0','60','120','180','240','300','360'});
    xlabel('degrees'); ylabel('trial');
    subplot(2,2,2);
    imagesc(stack(logical(idx),:));
    colorbar;
    title('Reward');
    set(gca,'xtick',0:size(stack,2)/6:size(stack,2));
    set(gca,'xticklabel',{'0','60','120','180','240','300','360'});
    xlabel('degrees'); ylabel('trial');

% %     theta=2*pi:-2*pi/size(stack,1):0;
%     theta=deg2rad(unit_rot(trials_ts(2):-1:trials_ts(1)));
%     count=1;
%     while(length(theta)~=size(stack,2))
%         theta=deg2rad(unit_rot(trials_ts(2+count):-1:trials_ts(1+count)));
%         count=count+1;
%     end
    subplot(2,2,3);
%     polarplot(theta,sum(stack(~idx,:))./size(stack(~idx,:),2)); hold on;
    polarplot(deg2rad(theta(1:end-1)),sum(stack(~idx,:))./size(stack(~idx,:),2)); hold on;
    polarplot([0 mean_vect_theta(1,n)],[0 mean_vect_rho(1,n)]); hold off;
    title(['p = ' mat2str(pval(n))]);
    subplot(2,2,4);
%     polarplot(theta,sum(stack(logical(idx),:))./size(stack(logical(idx),:),2)); hold on;
    polarplot(deg2rad(theta(1:end-1)),sum(stack(logical(idx),:))./size(stack(logical(idx),:),2)); hold on;
    polarplot([0 mean_vect_theta(1,n)],[0 mean_vect_rho(1,n)]); hold off;
    title(['p = ' mat2str(pval(n))]);
    
    supertitle({fn, ['n = ' mat2str(n)]});
    colormap hot;
    
    %PSTH
%     close all
%     vel=360/diff(frame_ts(trials_ts(1:2)));
%     err=std(stack(logical(idx),:))./sqrt(size(stack(logical(idx),:),1));
%     errorbar(([0:360/bin:360-360/bin]-270)./vel,flip(mean(stack(logical(idx),:)),2),err);
%     ylabel('average deconvolved')
%     xlabel('peristimulus time (s)')