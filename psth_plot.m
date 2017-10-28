function psth_plot(hObj,evt)
    global n
    global h
    global den
    global bins
    global trials
    global frame_ts
    global unit_pos
    
    switch evt.Key
        case 'rightarrow'
            n=n+1;
            if n>size(den,2)
                n=size(den,2);
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
    
    vel_thres=0.05;
    idx=diff(unit_pos);
    idx=find(idx>vel_thres | idx<-vel_thres);
    pos=unit_pos(idx);
    signal=den(idx,n);
    ft=frame_ts(idx);
        
    psth=zeros(length(trials)-1,bins);
    pre=find(ft>=trials(1),1);
    for i=2:length(trials)
        next=find(ft>=trials(i),1);
        p=pos(pre:next);
%         ratio=smooth(den(pre:next,n),31);
        ratio=signal(pre:next);
%         ratio=den(pre:next,n);
        idx=min(p):range(p)/bins:max(p);
        for j=1:bins
            psth(i-1,j)=mean(ratio(p>=idx(j) & p<=idx(j+1)));
        end
        pre=next;
    end
    psth=Smooth(psth,[0 2]);
    imagesc(psth);
    title(['n = ' mat2str(n)]);
    colormap hot
    ylabel('trials')
    xlabel('distance (arbitrary)')
    colorbar
%     hold on
%     plot([bins/4 bins/4],get(gca,'ylim'),'c','markersize',20);
%     plot([bins/2 bins/2],get(gca,'ylim'),'c','markersize',20);
%     plot([3*bins/4 3*bins/4],get(gca,'ylim'),'c','markersize',20);
%     idx=get(gca,'ylim');
%     idx=idx(1)+0.2;
%     text(5,idx,'light 1','color','cyan');
%     text(bins/4+5,idx,'tunnel 1','color','cyan');
%     text(bins/2+5,idx,'light 2','color','cyan');
%     text(3*bins/4+5,idx,'tunnel 2','color','cyan');