function pc_list=gmmtest(o_psth,o_stack,plotFlag)
pc_list=[];
count=1;
for k=1:size(o_psth,3)
    psth=o_psth(:,:,k);
    stack=o_stack(:,k);
    
    peaks=findpeaks(stack);
    peaks=peaks.loc(stack(peaks.loc)>(range(stack)*.3+min(stack)));
    nComp=length(peaks)+1; % determine the number of gaussian components by counting the number of peaks exceeding 30% of the FR range
    idx=psth>0;
    
    initial={};
    
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
    try
        fit=fitgmdist(x,nComp,'covariancetype','diagonal','sharedcovariance',false,'options',statset('maxiter',500));
    catch
        fit=fitgmdist(x,nComp,'covariancetype','diagonal','sharedcovariance',false,'options',statset('maxiter',500),'regularize',0.5);
    end
    thres=sqrt(chi2inv(0.95,2)); % cluster at 95% CI
    mahalDist=mahal(fit,x);
    
    clusters=double(mahalDist<=thres);
    for i=1:size(clusters,2)
        clusters(:,i)=clusters(:,i).*i;
    end
    clusters=sum(clusters,2);

    [~,noise]=min(fit.mu(:,2)); % find noise cluster
    idx=ones(1,fit.NumComponents);
    idx(noise)=0;

    sd=sqrt(fit.Sigma(:,1,:)); % s.d. over x
    sd=reshape(sd,1,size(sd,3));
    sd=4.*sd<length(stack); % components with +/-2*s.d. smaller than the width of the track

    idx=find(~(idx.*sd)); % bad clusters

%     clusters=cluster(fit,x);

    for i=idx % get rid of bad clusters
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
        x1 = linspace(min(x(:,1)) - 2,max(x(:,1)) + 2,500);
        x2 = linspace(min(x(:,2)) - 2,max(x(:,2)) + 2,500);
        [x1grid,x2grid] = meshgrid(x1,x2);
        X0 = [x1grid(:) x2grid(:)];
        mahalDist=mahal(fit,X0);
        if count>25
            figure;
            count=1;
        end
        subplot(5,5,count);
        hold on
        for i=1:nComp
            idx=mahalDist(:,i)<=thres;
            plot(X0(idx,1),X0(idx,2),'.','Color','c','MarkerSize',1);
        end
        gscatter(x(:,1),x(:,2));
        xlim([0 length(stack)]);
        title(num2str(passed));
        count=count+1;
    end
end