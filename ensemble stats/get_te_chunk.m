function te=get_te_chunk(deconv,order,precision)
% x -> y embedding

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

edges=[-1e6 min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
edges(end)=1e6;
deconv(isnan(deconv))=min(min(deconv))-1;
te=zeros(size(deconv,2));

dq=parallel.pool.DataQueue;
f=waitbar(0,'performing high-dimensional transfer entropy...');
afterEach(dq,@updateBar);
N=size(deconv,2);
p=1;

chunked=zeros(size(deconv,1)-order,size(deconv,2)*(1+order)); %allocate orders as chunks into memory for faster access
for i=0:order
    chunked(:,1+i*size(deconv,2):(i+1)*size(deconv,2))=deconv(order+1-i:end-i,:);
end

% chunked=gpuArray(chunked);
% te=gpuArray(te);

chunked=discretize(chunked,edges);
l=size(chunked,1);
s=size(deconv,2);
for i=1:s
    te(i,:)=arrayfun(@(x) embed_dimensions_gpu(i,x,s,l,order,chunked), 1:s);
    send(dq,i);
end

% te=gather(te);

close(f);

    function updateBar(~)
        waitbar(p/N, f);
        p=p+1;
    end

end
