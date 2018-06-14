function deconv=chain_align(deconv,D,lags,assemblies,id)
% realign deconv for assembly 'id'

deconv=deconv(:,assemblies{id});
D=D(assemblies{id},assemblies{id},:);

[D,idx]=min(D,[],3);

% D(diag(true(1,length(D))))=1;

% lags=reshape(lags,1,1,length(lags));
% lags=repmat(lags,size(D,1),size(D,2),1);
lags=lags(idx);

D(diag(true(1,length(D))))=1;

[~,chain]=min(min(D));
D(chain,:)=1;
lags_chain=[];
while length(chain)<length(D)
    [~,next]=min(D(:,chain(end)));
    D(next,:)=1;
    lags_chain=[lags_chain lags(next,chain(end))];
    chain=[chain next];
end

lags_chain=cumsum(lags_chain);
r=range(lags_chain);
shift=ceil(r/2)+min(lags_chain);
lags_chain=lags_chain-shift;

new_deconv=zeros(size(deconv,1)-r,size(deconv,2));

new_deconv(:,chain(1))=deconv(ceil(r/2)-shift+1:end-floor(r/2)-shift,chain(1));

for i=1:length(lags_chain)
    new_deconv(:,chain(i+1))=deconv(ceil(r/2)-lags_chain(i)+1:end-floor(r/2)-lags_chain(i),chain(i+1));
end

deconv=new_deconv;