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
% colormap gray
colormap(get_colour('black'));
hold on

l=cell(1,length(assemblies));
for i=1:length(assemblies)
    temp=deconv;
    temp(:,setxor(assemblies{i},1:size(deconv,2)))=0;
    if iscell(sce)
        temp(:,assemblies{i})=temp(:,assemblies{i}).*sce{i};
    else
        temp(~sce(:,i),:)=0;
    end
    temp=temp(:,ordered);
    [x,y]=find(temp);
    plot(x,y,'.');
    l{i}=['Ensemble ' num2str(i)];
end

legend(l,'location','southeast');
xlabel('frame')
ylabel('ordered neuron no.')