function monte_carlo_assemblies(deconv)
% Assemblies detection algorithm using methods in independent study

deconv=ca_filt(deconv); % get rid of "noise"
r=deconv>0; % binarize spiking probability matrix

pn=sum(r)./size(r,1);
pn=repmat(pn,size(r,1),1);
% P=r.*pn+(~r).*(1-pn);
% P=prod(P,2);

% new monte carlo simulation
number_sims=1000;
chunk=1000; % chunks to process per iteration
P=zeros(1,number_sims*size(r,1));
C=P;
pn=repmat(pn,chunk,1);
for n=1:number_sims/chunk
    temp=r;
    
    temp=repmat(temp,1,chunk);
    
    shift=randi(size(temp,1),1,size(temp,2));
    temp=mat_circshift(double(temp),shift,1);
    
%     for i=1:size(r,2)
%         temp(:,i)=circshift(temp(:,i),randi(size(temp,1)));
%     end

    temp=reshape(temp,[],size(r,2));

    P((n-1)*size(temp,1)+1:n*size(temp,1))=prod(temp.*pn+(~temp).*(1-pn),2);
    C((n-1)*size(temp,1)+1:n*size(temp,1))=sum(temp,2);
end