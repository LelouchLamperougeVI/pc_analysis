function cell_plot(hObj,evt)
    global n
    global h
    global trials_ts
    global deconv
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
        
    stack=zeros(max(diff(trials_ts))+1,length(trials_ts));
    for i=1:length(trials_ts)-1
        try
            stack(:,i)=smooth(deconv(trials_ts(i):trials_ts(i+1),n),21);
        catch
            stack(1:end-1,i)=smooth(deconv(trials_ts(i):trials_ts(i+1),n),21);
        end
    end

    subplot(2,2,1);
    imagesc(stack(:,1:2:end)');
    subplot(2,2,2);
    imagesc(stack(:,2:2:end)');

    theta=0:2*pi/size(stack,1):2*pi;
    subplot(2,2,3);
    polarplot(theta(1:end-1),sum(stack(:,1:2:end),2)./size(stack,1));
    subplot(2,2,4);
    polarplot(theta(1:end-1),sum(stack(:,2:2:end),2)./size(stack,1));
    
    title(['n = ' mat2str(n)]);