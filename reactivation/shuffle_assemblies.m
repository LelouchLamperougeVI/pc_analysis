function sd_series=shuffle_assemblies(deconv,plotFlag)
% Traditional method used by Yuste and Cossart
% Same principal as monte_carlo_assemblies

if nargin<2
    plotFlag=false;
end

deconv=ca_filt(deconv); % get rid of "noise"
r=deconv>0; % binarize spiking probability matrix

number_sims=1000;
C=zeros(1,number_sims*size(r,1));
for i=1:number_sims
    temp=r;
    
    shift=randi(size(temp,1),1,size(temp,2));
    temp=mat_circshift(double(temp),shift);

    C((i-1)*size(temp,1)+1:i*size(temp,1))=sum(temp,2);
end

[mu,sigma]=normfit(C);
sd_series=(sum(r,2)-mu)./sigma;
% sd_series(isnan(sd_series))=4;

if plotFlag
    figure;
    ax1=subplot(2,1,1);
    mat_scatter(r,'k');
    ax2=subplot(2,1,2);
    plot(sd_series);
    linkaxes([ax1,ax2],'x');
    ylabel('s.d.'); xlabel('frame');
end