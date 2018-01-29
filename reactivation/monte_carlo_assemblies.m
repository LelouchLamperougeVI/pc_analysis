function sd_series=monte_carlo_assemblies(deconv, plotFlag)
% Assemblies detection algorithm using methods in independent study
% Returns number of standard deviations from the mean of each cell-count distribution

if nargin<2
    plotFlag=false;
end

deconv=ca_filt(deconv); % get rid of "noise"
r=deconv>0; % binarize spiking probability matrix

pn=sum(r)./size(r,1);
pn=repmat(pn,size(r,1),1);

% new monte carlo simulation
number_sims=1000;
P=zeros(1,number_sims*size(r,1));
C=P;
for i=1:number_sims
    temp=r;
    
    shift=randi(size(temp,1),1,size(temp,2));
    temp=mat_circshift(double(temp),shift);

    P((i-1)*size(temp,1)+1:i*size(temp,1))=prod(temp.*pn+(~temp).*(1-pn),2);
    C((i-1)*size(temp,1)+1:i*size(temp,1))=sum(temp,2);
end

means=zeros(2,size(r,2)+1);
for i=1:size(r,2)+1
    means(:,i)=lognfit(P(C==i-1));
end

P=prod(r.*pn+(~r).*(1-pn),2);
C=sum(r,2);
sd_series=(log(P)'-means(1,C+1))./means(2,C+1);
sd_series(isnan(sd_series))=4;

if plotFlag
    figure;
    ax1=subplot(2,1,1);
    mat_scatter(r,'k');
    ax2=subplot(2,1,2);
    plot(sd_series);
    linkaxes([ax1,ax2],'x');
    ylabel('s.d.'); xlabel('frame');
end