function [assemblies,R]=lopes_pca(deconv,smooth,plotFlag)
% PCA cell-assemblies detection methods adapted from Lopes-dos-Santos et
% al. 2011
% Inputs:
%   deconv: raw data with R samples and C neurons
%   smooth: smoothing factor
%   plotFlag: set true for relevant figures and verbose
% Outputs:
%   assemblies: cell array of assembly neurons
%   R: activation strength for R samples and C assemblies

if nargin<3
    plotFlag=false;
end

deconv=fast_smooth(deconv,smooth);
deconv=(deconv-mean(deconv))./std(deconv); % zscore each neuron's spike activity
corrM=corr(deconv); % compute autocorrelation matrix
[V,e]=eig(corrM); % get eigenvalues and vectors/PCs
e=sum(e); % vectorize eigenvalues matrix

[thres,epdf]=MPdist(size(deconv,1),size(deconv,2),min(e):0.01:max(e)); % get upper/lower thresholds for eigenvalue significance
if plotFlag
    figure;
    histogram(e,50,'normalization','pdf');
    hold on;
    plot(min(e):0.01:max(e),real(epdf),'r');
    xlabel('eigenvalue');
    ylabel('frequency');
end

assemblies=e>thres(2); % identify PC's associated with assemblies
numNeurons=sum(e<thres(1) | e>thres(2)); % determine number of neurons participating in assemblies
assembly_space=V(:,assemblies); % PCs with significant eigenvalues
interaction=assembly_space*assembly_space'; % interaction matrix
vect_lengths=arrayfun(@(x) norm(assembly_space(x,:)),1:size(assembly_space,1));
[~,idx]=sort(vect_lengths);
assembly_neurons=idx(end-numNeurons+1:end); % neurons with biggest norms in assembly space
interaction=interaction(assembly_neurons,:); % remove non-assembly neurons from interaction matrix
interaction=interaction(:,assembly_neurons);
diagonal=repmat(vect_lengths(assembly_neurons)',1,size(interaction,1));
interaction=interaction./diagonal'; % final interaction matrix is the scalar projection of neurons onto each other in assembly space

linearized_interaction=reshape(interaction,[],1);
idx=kmeans(linearized_interaction,2); % k-means clustering to set threshold for significant interactions
idx=idx-1;
idx=logical(idx);
if mean(linearized_interaction(idx)) < mean(linearized_interaction(~idx))
    idx=~idx;
end

if plotFlag
    figure;
    hist(linearized_interaction,1000);
    hold on;
    plot([min(linearized_interaction(logical(idx))) min(linearized_interaction(logical(idx)))],[0 max(get(gca,'ylim'))],'r--');
    title('Interaction Threshold');
    xlabel('projection');
    ylabel('counts');
end

identicals=idx==idx'; % remove identical rows for BIM
for i=1:size(interaction,1)
    for j=1:size(interaction,2)
        ident(i,j)=prod(diag(identicals(size(interaction,1)*(i-1)+1:size(interaction,1)*i,size(interaction,1)*(j-1)+1:size(interaction,1)*j)));
    end
end
diagonal=zeros(size(ident,1));
for i=1:size(ident,1)-1
    diagonal=diagonal+diag(ones(1,size(ident,1)-i),i);
end
ident=ident.*diagonal;
ident=logical(sum(ident,2));

BIM=reshape(idx,[],size(interaction,1)); % generate binary interaction matrix
BIM(ident,:)=[]; % remove identical rows

count=sum(BIM,2);
[~,idx]=sort(count);
BIM=BIM(idx,:); % reorder BIM rows by increasing number of interactions

BIM_index=ones(size(interaction,1),1); % rows index for BIM
BIM_index(ident)=0;
BIM_index=find(BIM_index);
BIM_index=BIM_index(idx);

ALM=cell(size(interaction,1)); % clustering algorithm
L=0;
r=1;
while(r<=size(BIM,1))
    rr=BIM_index(r);
    if isempty(cell2mat(ALM(rr,:)))
        L=L+1;
        for i=find(BIM(r,:))
            for j=find(BIM(r,:))
                ALM(i,j)={[ALM{i,j} L]};
            end
        end
    end
%     ALM(BIM(r,:),BIM(r,:))={[cell2mat(ALM(BIM(r,:),BIM(r,:))) L]};
    r=r+1;
end

assemblies=cell(1,L);
for i=1:size(ALM,1)
    for j=cell2mat(ALM(i,i))
        assemblies(j)={sort([cell2mat(assemblies(j)) assembly_neurons(i)])};
    end
end

if plotFlag
    fprintf('Detected assemblies: \n');
    arrayfun(@(x) fprintf(['\t Assembly ' num2str(x) ':\t' mat2str(assemblies{x}) '\n']),1:length(assemblies));
end

R=zeros(size(deconv,1),length(assemblies)); % the math in the original paper was a bit dodgy (or maybe I'm just too dumb to understand) so I changed it a bit without affecting the outcomes
for i=1:length(assemblies)
    a=sum(assembly_space(assemblies{i},:))./norm(sum(assembly_space(assemblies{i},:)));
%     alpha=sum(repmat(a*assembly_space(assemblies{i},:),length(assemblies{i}),1).*assembly_space(assemblies{i},:),2);
    alpha=sum(repmat(a,size(assembly_space,1),1).*assembly_space,2);
    psi=alpha*alpha';
%     Z=deconv(:,assemblies{i});
    R(:,i)=sum(deconv*psi.*deconv,2);
end

if plotFlag
    figure;
    ax1=subplot(2,1,1);
    imagesc(deconv');
    colormap gray;
    ax2=subplot(2,1,2);
    plot(R);
    linkaxes([ax1,ax2],'x');
    labels=cell(1,size(R,2));
    for i=1:size(R,2); labels{i}=['Ensemble ' num2str(i)]; end
    legend(labels);
end
