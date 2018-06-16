function te=get_te_chunk(deconv,order,precision)
% x -> y embedding

deconv=ca_filt(deconv);
deconv=log(deconv);
deconv(isinf(deconv))=nan;
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');

edges=[-1e6 min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
edges(end)=1e6;
deconv(isnan(deconv))=min(min(deconv))-1;

deconv=single(deconv);

% te=gpuArray(zeros(size(deconv,2)));
te=zeros(size(deconv,2));

% dq=parallel.pool.DataQueue;
f=waitbar(0,'performing high-dimensional transfer entropy...');
% afterEach(dq,@updateBar);
% N=size(deconv,2)^2;
N=size(deconv,2);
p=1;

chunked=zeros(size(deconv,1)-order,size(deconv,2)*(1+order)); %allocate orders as chunks into memory for faster access
for i=0:order
    chunked(:,1+i*size(deconv,2):(i+1)*size(deconv,2))=deconv(order+1-i:end-i,:);
end

% % chunked=gpuArray(chunked);

chunked=discretize(chunked,edges);
l=size(chunked,1);
s=size(deconv,2);
% tmp=gpuArray(zeros(1,size(deconv,2)));
% tmp=zeros(1,size(deconv,2));
for i=1:s
%     tmp=zeros(1,size(deconv,2));
%     for j=1:s
%         tmp(j)=embed_dimensions_gpu(i,j,s,l,order,chunked);
%         send(dq,i);
%         updateBar;
%     end
    te(i,:)=arrayfun(@(x) embed_dimensions_gpu(i,x,s,l,order,chunked), 1:s);
%     te(i,:)=tmp;
    updateBar;
end

te=gather(te);

close(f);

%     function te=embed_dimensions(x,y)
%         % x -> y embedding
%         
%         idx=[y+s.*(0:order) x+s.*(1:order)];
%         p_joint=accumarray(chunked(:,idx),1)./l;
%         
%         idx=[y+s.*(1:order) x+s.*(1:order)];
%         p_cond_xy=accumarray(chunked(:,idx),1)./l;
%         p_cond_xy=permute(p_cond_xy,[ndims(p_cond_xy)+1 1:ndims(p_cond_xy)]); %inject singletons
%         p_cond_xy=p_joint./p_cond_xy;
%         
%         idx=y+s.*(0:order);
%         p_joint_yy=accumarray(chunked(:,idx),1)./l;
%         
%         idx=y+s.*(1:order);
%         p_joint_y=accumarray(chunked(:,idx),1)./l;
%         p_joint_y=permute(p_joint_y,[ndims(p_joint_y)+1 1:ndims(p_joint_y)]);
%         p_cond_y=p_joint_yy./p_joint_y;
%         
%         te=p_joint.*log2(p_cond_xy./p_cond_y);
%         te=sum(te(:),'omitnan');
%     end

    function updateBar(~)
        waitbar(p/N, f);
        p=p+1;
    end

end
