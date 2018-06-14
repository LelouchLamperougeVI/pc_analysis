function te=get_te(deconv,order,precision)
% x -> y embedding

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

edges=[-inf min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
edges(end)=inf;
deconv(isnan(deconv))=min(min(deconv))-1;

te=zeros(size(deconv,2));

dq=parallel.pool.DataQueue;
f=waitbar(0,'performing high-dimensional transfer entropy...');
afterEach(dq,@updateBar);
N=size(deconv,2);
p=1;

parfor i=1:size(deconv,2)
    tmp=zeros(1,size(deconv,2));
    for j=1:size(deconv,2)
        x=deconv(:,i);
        y=deconv(:,j);
        tmp(j)=embed_dimensions(x,y,edges,order);
    end
    te(i,:)=tmp;
    send(dq,i);
end

close(f);

    function updateBar(~)
        waitbar(p/N, f);
        p=p+1;
    end

end
