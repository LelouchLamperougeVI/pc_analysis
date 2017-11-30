function pc_list=gmmtest(o_psth,o_stack,plotFlag)
pc_list=[];
for k=1:size(o_psth,3)
    psth=o_psth(:,:,k);
    stack=o_stack(:,k);
    
    peaks=findpeaks(stack);
    peaks=peaks.loc(stack(peaks.loc)>(range(stack)*.3+min(stack)));
    nComp=length(peaks)+1; % determine the number of gaussian components by counting the number of peaks exceeding 30% of the FR range
    idx=psth>0;
    % thres=range(range(psth))*.05;
    % idx=psth>thres;
    x=[];
    trials=[];
    for i=1:size(psth,1)
        for j=1:size(psth,2)
            if idx(i,j)
                x=[x;j log(psth(i,j))];
    %             x=[x;j psth(i,j)];
                trials=[trials i];
            end
        end
    end
    fit=fitgmdist(x,nComp);

    [~,noise]=min(fit.mu(:,2)); % find noise cluster
    idx=ones(1,fit.NumComponents);
    idx(noise)=0;

    sd=sqrt(fit.Sigma(1,1,:)); % s.d. over x
    sd=reshape(sd,1,size(sd,3));
    sd=2.*sd<length(stack); % components with s.d. smaller than the width of the track

    idx=find(~(idx.*sd));

    clusters=cluster(fit,x);

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
end
    
    
if nargin>2 && plotFlag
    figure;
    h = gscatter(x(:,1),x(:,2));
    hold on; 
    fcontour(@(x1,x2)pdf(fit,[x1 x2]),[0 50 -18 3])
end