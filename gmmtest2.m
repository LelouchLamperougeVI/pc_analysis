function pc_list=gmmtest2(o_psth,o_stack,plotFlag)
pc_list=[];
count=1;
for k=1:size(o_psth,3)
    psth=o_psth(:,:,k);
    stack=o_stack(:,k);
    
    peaks=findpeaks(stack);
    peaks=peaks.loc(stack(peaks.loc)>(range(stack)*.3+min(stack)));
    nComp=length(peaks); % determine the number of gaussian components by counting the number of peaks exceeding 30% of the FR range
    
    idx=psth>0;
    x=[];
    trials=[];
    for i=1:size(psth,1)
        for j=1:size(psth,2)
            if idx(i,j)
                x=[x;j log(psth(i,j))];
                trials=[trials i];
            end
        end
    end
    % isolate noise region
    fit=fitgmdist(x(:,2),2);
    [~,idx]=min(fit.mu);
    thres=sqrt(chi2inv(0.95,2)); % cluster at 95% CI
    mahalDist=mahal(fit,x(:,2));
    noise=mahalDist(:,idx)<=thres;
    
    try
        fit=fitgmdist(x(~noise,1),nComp);
    catch
        fit=fitgmdist(x(~noise,1),nComp,'regularize',1);
    end
    sigma=reshape(sqrt(fit.Sigma).*4,size(fit.Sigma,3),1)>length(stack);
    mahalDist=mahal(fit,x(:,1));
    clusters=cluster(fit,x(:,1));
    
    clusters=logical(sum(mahalDist<=thres,2)).*~noise.*clusters;
    for i=find(sigma)' % get rid of bad clusters
        clusters(clusters==i)=0;
    end
    
    passed=0;
    for i=unique(clusters(clusters>0))'
        idx=clusters==i;
        idx=idx'.*trials;
        num_trials=length(unique(trials(logical(idx))));
        if num_trials>size(psth,1)/3
            passed=1;
        end
    end
    
    if passed
        pc_list=[pc_list k];
    end
    
    if nargin>2 && plotFlag
        if count>25
            figure;
            count=1;
        end
        subplot(5,5,count);
        hold on
        plot(0:length(stack)-1,stack);
        plot(0:0.1:length(stack),pdf(fit,(0:0.1:length(stack))'));
%         xlim([0 length(stack)]);
        title(num2str(passed));
        count=count+1;
    end
end