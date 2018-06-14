function te=embed_dimensions(x,y,edges,order)
% x -> y embedding

dims=zeros(length(y)-order,2*order+1);
dims(:,1)=y(order+1:end);
for i=1:order
    dims(:,i*2)=y(order+1-i:end-i);
    dims(:,i*2+1)=x(order+1-i:end-i);
end
dims=discretize(dims,edges); %can't believe I didn't discover this shit sooner

p_joint=accumarray(dims,1)./size(dims,1);

p_joint_xy=accumarray(dims(:,2:end),1)./size(dims,1);
p_joint_xy=permute(p_joint_xy,[ndims(p_joint_xy)+1 1:ndims(p_joint_xy)]); %inject singletons
p_cond_xy=p_joint./p_joint_xy;

p_joint_yy=accumarray(dims(:,[1 2:2:end]),1)./size(dims,1);
idx=zeros(1,order*2);
idx(1:2:end)=2:ndims(p_joint_yy);
idx(2:2:end)=ndims(p_joint_yy)+(1:order);
p_joint_yy=permute(p_joint_yy,[1 idx]);

p_joint_y=accumarray(dims(:,2:2:end),1)./size(dims,1);
idx=zeros(1,order*2);
idx(1:2:end)=1:ndims(p_joint_y);
idx(2:2:end)=ndims(p_joint_y)+(2:order+1);
p_joint_y=permute(p_joint_y,[ndims(p_joint_y)+1 idx]);
p_cond_y=p_joint_yy./p_joint_y;

te=p_joint.*log2(p_cond_xy./p_cond_y);
% te(isinf(te))=0;
te=sum(te(:),'omitnan');