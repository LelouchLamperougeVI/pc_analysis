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
