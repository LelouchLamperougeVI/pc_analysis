function te=get_te(deconv,order,precision)
% x -> y embedding

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

edges=[-inf min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
edges(end)=inf;
deconv(isnan(deconv))=min(min(deconv))-1;

deconv=single(deconv);

te=zeros(size(deconv,2));

dq=parallel.pool.DataQueue;
f=waitbar(0,'performing high-dimensional transfer entropy...');
afterEach(dq,@updateBar);
N=size(deconv,2)^2;
p=1;

for i=1:size(deconv,2)
    tmp=zeros(1,size(deconv,2));
    x=deconv(:,i);
    parfor j=1:size(deconv,2)
        y=deconv(:,j);
        tmp(j)=embed_dimensions_sparse(x,y,edges,order);
        send(dq,i);
    end
    te(i,:)=tmp;
end

close(f);

    function updateBar(~)
        waitbar(p/N, f);
        p=p+1;
    end

end
