function cossart_sequences(deconv,tcs,smooth)
% Sequence detection from calcium signals adapted from Malvache et al. 2016

deconv=ca_filt(deconv); % not sure if this step is even necessary
[assemblies,R]=lopes_pca(deconv,smooth); % cell-assemblies detection with Lopes-dos-Santos algorithm

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
    
    raw=double(tcs.ratio(:,assembly));
    for k=1:size(events,2) % get sequence for each event
        w=raw(events(1,k):events(2,k),:);
        template=median(w);
        
    end
end