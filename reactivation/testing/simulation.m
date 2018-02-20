%% percent accuracy as function of network size and time samples

m=520;
n=520;

A=zeros(20,20);

for i=500:m
  for j=500:n
%     q=bin_poisson_sim(i,j,'R',0.05); %estimated from real data
    q=bin_poisson_sim(i,j,'R',prob(1:j)); %estimated from real data
    assemblies=bin_pca(q);
    A(i-499,j-499)=assemblies.accuracy;
  end
  i
end


%% detection accuracy as a function of firing rate and activation probability

m=8000;
n=150;
numAssemblies=10;
neuAssembly=10;
assemblies=zeros(numAssemblies,n);
for i=1:numAssemblies
  assemblies(i,neuAssembly*(i-1)+1:neuAssembly*i)=1;
end
R=0.001:0.001:0.2;
probability=0.001:0.001:0.05;

A=zeros(length(R),length(probability));
parfor i=1:length(R)
  temp=zeros(1,length(probability));
  for j=1:length(probability)
    q=bin_poisson_sim(m,n,'R',R(i),'prob',probability(j),'assemblies',assemblies);
    a=bin_pca(q);
    temp(j)=a.numAssemblies/numAssemblies;
  end
  A(i,:)=temp;
end

%% detection count accuracy as a function of assembly size

m=8000;
n=150;
numAssemblies=1;
neuAssembly=2:100;
R=0.05;
probability=0.001:0.001:0.05;

A=zeros(length(neuAssembly),length(probability));
parfor i=1:length(neuAssembly)
  temp=zeros(1,length(probability));
  for j=1:length(probability)
    assemblies=zeros(numAssemblies,n);
    assemblies(i,1:neuAssembly(i))=1;
    q=bin_poisson_sim(m,n,'R',R,'prob',probability(j),'assemblies',assemblies);
    a=bin_pca(q);
    temp(j)=a.numAssembliesNeurons/neuAssembly(i);
  end
  A(i,:)=temp;
end

%% accuracy as a function of background firing rate and missing members

m=8000;
n=150;
numAssemblies=1;
neuAssembly=50;
R=0.01:0.01:0.5;
probability=0.05;
miss=0:25;

A=zeros(length(R),length(miss));
    assemblies=zeros(numAssemblies,n);
    assemblies(1,1:neuAssembly)=1;
parfor i=1:length(R)
  temp=zeros(1,length(probability));
  for j=1:length(miss)
    q=bin_poisson_sim(m,n,'R',R(i),'prob',probability,'assemblies',assemblies,'miss',miss(j));
    a=bin_pca(q);
    temp(j)=a.numAssembliesNeurons/neuAssembly;
  end
  A(i,:)=temp;
end

%% sequences
run_epochs=noRun(behavior.speed_raw);
run_epochs=behavior.speed_raw>run_epochs;
jitter=4;
run_epochs=logical(conv(run_epochs,ones(1,jitter),'same'));
deconv(run_epochs,:)=[];
ratio_model(run_epochs,:)=[];

smooth=3;
q=logical(ca_filt(deconv));
q=adjacency_matrix(q,smooth);
% q=bin_q(q,10);
[assemblies,R]=lopes_pca(q,0,1);

perc=95; %percentile
idx=prctile(R,perc);
idx=R>=idx;
idx=get_head(idx);
% idx=get_middle(idx);

interp_size=10;
t=1:1/interp_size:size(q,1);
interp_ts=zeros(length(t),size(q,2));
for i=1:size(q,2)
    interp_ts(:,i)=interp1(1:size(q,1),ratio_model(:,i),t,'pchip');
end
interp_ts=detrend(interp_ts);
interp_ts=zscore(interp_ts);

sequence_mat=zeros(size(interp_ts));
e=logical(sum(idx,2));
e=get_head(e);
e=find(e);
idx=adjacency_matrix(idx,smooth);
idx=logical(idx);
for i=1:length(e)
    assembly_cells=idx(e(i),:);
    assembly_cells=unique(cell2mat(assemblies(assembly_cells)));
    
    [sequence,no_spike]=get_sequence(assembly_cells,interp_ts,e(i),t,smooth*2);
    
    sequence(no_spike)=[];
    assembly_cells(no_spike)=[];
    sequence_mat(sequence,assembly_cells)=1;
end

%%
R(~idx)=0;
[~,e]=max(R');
e(~logical(sum(idx,2)))=0;

sequences=cossart_sequences2(raw,e,assemblies);

for seq=1:8
    thres=200; %frames
    lol=e(e>0);
    lol=find(lol==seq);
    lol([false diff(find(e==seq))<thres])=[];
    lol=sequences(lol,:);
    match=lol;
    [~,order]=sort(mean(lol));
    lol=lol(:,order);
    lol=lol+repmat(cumsum(max(lol')'),1,size(lol,2));
    temp=size(lol,1);
    lol=reshape(lol',1,[]);
%     plot(lol,repmat(1:size(sequences,2),1,temp),'.');

    sampling_rate=19.1; %fps
    lag_window=0.5; %seconds
    e_window=round(lag_window*sampling_rate/2);

    stack=analysis.raw_stack;
    [~,template]=max(stack);
    template=template./50.*mean(diff(behavior.trials));

    match=mean(match)./(2*e_window+1).*lag_window;
    template(isnan(match))=[];
    match(isnan(match))=[];

    for cf=1:20
        coef(cf)=corr(match',template'./cf);
    end
    hold on;
    plot(coef);
end

[~,lol]=sort(mean(sequences));
lol=sequences(:,lol);
lol=lol+repmat(cumsum(max(lol')'),1,size(lol,2));
colors='bkrgycmw';
idx=e(e>0);
hold on
for i=1:size(lol,1)
    c=colors(idx(i));
    plot(lol(i,:),1:size(sequences,2),['.' c]);
end







