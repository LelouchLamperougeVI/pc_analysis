function sequences=cossart_sequences(deconv,tcs,smooth)
% Sequence detection from calcium signals adapted from Malvache et al. 2016

deconv=ca_filt(deconv); % not sure if this step is even necessary
[assemblies,R]=lopes_pca(deconv,smooth); % cell-assemblies detection with Lopes-dos-Santos algorithm

sequences=cell(1,length(assemblies));
for i=1:length(assemblies)
    assembly=assemblies{i};
    r=R(:,i); % FYI this distribution appears lognormal
    param=lognfit(r); % fit activation strength to lognormal distribution
    idx=r>(2*param(2)+param(1)); % reactivation epochs defined as those with 2 s.d. above the mean reactivation strength
    events=zeros(1,length(idx)); % find event onset/offset times
    if idx(1)==1
        events(1)=1;
    end
    if idx(end)==1
        events(end)=-1;
    end
    idx=diff(idx);
    events(find(idx==1)+1)=events(find(idx==1)+1)+1;
    events(idx==-1)=events(idx==-1)-1;
    events=[find(events==1);find(events==-1)];
%     events=events+[-5;5]; % give it a +/- 5 samples wiggle room
    events(events<1)=1;
    events(events>length(r))=length(r);
    
    raw=double(tcs.ratio(:,assembly));
    sequence=zeros(length(assembly),size(events,2));
    for k=1:size(events,2) % get sequence for each event
        w=raw(events(1,k):events(2,k),:);
        template=median(w,2);
        vertex=zeros(1,size(w,2));
        for b=1:size(w,2)
            [covariance,lags]=xcov(w(:,b),template);
            p=polyfit(lags,covariance',2);
            vertex(b)=-p(2)/2/p(1);
            if max(covariance)<0.6
                vertex(b)=NaN;
            end
        end
        [~,idx]=sort(vertex);
        idx(isnan(vertex))=[];
        sequence(idx,k)=1:length(idx);
    end
    sequences{i}=sequence;
end