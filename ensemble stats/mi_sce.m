function [sce,loss]=mi_sce(assemblies,deconv,precision,max_it,jitter)
if nargin<4
    max_it=50;
end
if nargin<5
    jitter=50;
end

% deconv=logical(deconv);


sce=zeros(size(deconv,1),length(assemblies));
% idx=zeros(size(deconv,1),length(assemblies));
loss=cell(1,length(assemblies));
for i=1:length(assemblies)
    temp=deconv(:,assemblies{i});
%     shuffled_d=mat_circshift(temp,randi(size(temp,1),1,size(temp,2)));
%     [~,cutoff]=get_mi(shuffled_d);
    [~,cutoff]=get_mi(temp,precision);
    cutoff=triu(cutoff,1);
    cutoff=mean(cutoff(cutoff~=0));
    
    H=get_entropy(deconv,precision);
    [~,idx]=sort(H,'descend');
    temp=temp(idx,:);
    
    cost=0;
    converge=false;
    count=1;
%     while count<=size(temp,1)
%     while ~converge && count<=size(temp,1)
    while cost(end)<cutoff && ~converge && count<=size(temp,1)
        try
            [~,d]=get_mi(temp(1:count,:),precision);
            d=triu(d,1);
            d=mean(d(d~=0));
            d(isnan(d))=-inf;
            cost=[cost d];
%         try
            if cost(end)<=mean(cost(end-max_it))
                converge=true;
            end
        catch
        end
        count=count+1;
    end
    cost(1)=[];
    loss{i}=cost;
    
    sce(idx(1:count-1),i)=1;
    idx=find(sce(:,i));
    idx2=find(diff(idx)<=jitter);
    for j=1:length(idx2)
        sce(idx(idx2(j)):idx(idx2(j)+1),i)=1;
    end
    
    figure;
    ax1=subplot(2,3,1:2);
    imagesc(deconv(:,assemblies{i})');
    ylabel('neuron no.')
    ax2=subplot(2,3,4:5);
    plot(sce(:,i));
    xlabel('frame')
    ylim([-1 2])
    yticks([0 1]);
    yticklabels({'idle','activation'});
    linkaxes([ax1 ax2],'x');
    subplot(2,3,[3 6]);
    plot([loss{i}]);
    hold on
    plot(1:length([loss{i}]),ones(1,length([loss{i}])).*cutoff);
    xlabel('iteration')
    ylabel('d')
    legend({'d','complete sample d'},'location','southeast');
end