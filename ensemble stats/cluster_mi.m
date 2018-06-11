function [assemblies,Z]=cluster_mi(varargin)
% assemblies membership assignment using agglomerative clustering of
% MI distance metric
% Inputs:
%   'method':       time-delayed MI ('tdmi') or simply MI ('mi' - default)
%   'sig':          significance threshold (default 5%)
%   'prune':        prune small ensembles; minimum # members
%   'shuffle':      number of shuffles (default 500)
%   'precision':    number of bins for FR (default 5; not required for binary matrix)
%   'max_delay':    maximum delay between ensemble events (for tdmi only)(default 5)
%   'plotFlag':     plot relevant figures

deconv=varargin{1};
[method,max_delay,sig,prune,shuffle,precision,plotFlag]=parse_input(varargin);
if islogical(deconv)
    deconv=double(deconv);
end

if strcmpi(method,'mi')
    [I,D]=get_mi(deconv,precision);
else
    [I,D]=tdmi(deconv,precision,max_delay);
    D=min(D,[],3);
    I=max(I,[],3);
end

cutoff=[];
for i=1:shuffle
    shuffled_d=mat_circshift(deconv,randi(size(deconv,1),1,size(deconv,2)));
%     if strcmpi(method,'mi')
        [~,shuffled_d]=get_mi(shuffled_d,precision);
%     else
%         [~,shuffled_d]=tdmi(shuffled_d,precision,max_delay);
%         shuffled_d=min(shuffled_d,[],3);
%     end
    shuffled_d=triu(shuffled_d,1);
    cutoff=[cutoff;shuffled_d(~~shuffled_d)];
end
cutoff=prctile(cutoff,sig);

D=squareform(D);
Z=linkage(D,'average');
clusters=cluster(Z,'criterion','distance','cutoff',cutoff);
D=squareform(D);

idx=1:length(clusters);
% assemblies=cell(1,max(clusters));
count=1;
for i=1:max(clusters)
    assemblies{count}=idx(clusters==i);
    temp=D(assemblies{count},assemblies{count});
    temp=triu(temp,1);
    temp=temp(~~temp);
    if mean(temp)<cutoff && length(assemblies{count})>=prune
        count=count+1;
    else
        assemblies(count)=[];
    end
end


if plotFlag
    fprintf('Detected assemblies: \n');
    arrayfun(@(x) fprintf(['\t Assembly ' num2str(x) ':\t' mat2str(assemblies{x}) '\n']),1:length(assemblies));
    figure
    ax1=subplot(2,2,2);
    [~,order]=sort(clusters);
    dendrogram(Z,0,'colorthreshold',cutoff,'reorder',order);
    axis square
    ylabel('d')
    ax2=subplot(2,2,4);
    imagesc(I(order,order));
    axis square
    xlabel('neuron no.')
    ylabel('neuron no.')
    ax3=subplot(2,2,3);
    dendrogram(Z,0,'colorthreshold',cutoff,'reorder',order(end:-1:1),'orientation','left');
    axis square
    xlabel('d')
    linkaxes([ax1,ax2],'x');
    linkaxes([ax2,ax3],'y');
    
    %     figure;
    %     imagesc(deconv(:,cell2mat(assemblies))');
    %     colormap gray
end


function [method,max_delay,sig,prune,shuffle,precision,plotFlag]=parse_input(inputs)
precision=5;
plotFlag=false;
sig=5;
prune=0;
shuffle=500;
method=1;
max_delay=5;

idx=2;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 'prune'
            idx=idx+1;
            prune=inputs{idx};
        case 'shuffle'
            idx=idx+1;
            shuffle=inputs{idx};
        case 'precision'
            idx=idx+1;
            precision=inputs{idx};
        case 'method'
            idx=idx+1;
            method=inputs{idx};
        case 'max_delay'
            idx=idx+1;
            method=inputs{idx};
        case 'plotflag'
            idx=idx+1;
            plotFlag=inputs{idx};
        otherwise
    end
    idx=idx+1;
end

