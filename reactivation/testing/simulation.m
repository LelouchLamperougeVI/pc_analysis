%% percent accuracy as function of network size and time samples

m=520;
n=520;

A=zeros(20,20);

for i=500:m
  for j=500:n
    q=bin_poisson_sim(i,j,'R',0.05); %estimated from real data
    assemblies=bin_pca(q);
    A(i,j)=assemblies.accuracy;
  end
  i
end
