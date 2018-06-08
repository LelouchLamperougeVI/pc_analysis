function h=plot_assemblies(assemblies,sce,deconv,ordered)
if nargin<4
    ordered=cell2mat(assemblies);
    ordered=[ordered setxor(1:size(deconv,2),ordered)];
end

h=figure;
if islogical(deconv)
    imagesc(deconv(:,ordered)');
else
    imagesc(fast_smooth(deconv(:,ordered),5)');
end
colormap gray
hold on

l=cell(1,length(assemblies));
for i=1:length(assemblies)
    temp=deconv;
    temp(~sce(:,i),:)=0;
    temp(:,setxor(assemblies{i},1:size(deconv,2)))=0;
    temp=temp(:,ordered);
    [x,y]=find(temp);
    plot(x,y,'.');
    l{i}=['Ensemble ' num2str(i)];
end

legend(l,'location','southeast');
xlabel('frame')
ylabel('ordered neuron no.')