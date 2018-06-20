function te=embed_dimensions_gpu(x,y,s,l,order,chunked)
% x -> y embedding

idx=[y+s.*(0:order) x+s.*(1:order)];
p_joint=accumarray(chunked(:,idx),1)./l;

p_cond_xy=sum(p_joint);
p_cond_xy=p_joint./p_cond_xy;

for i=order+2:ndims(p_joint)
    p_joint_yy=sum(p_joint,i);
end

p_joint_y=sum(p_joint_yy,1);
p_cond_y=p_joint_yy./p_joint_y;

te=p_joint.*log2(p_cond_xy./p_cond_y);
te=sum(te(:),'omitnan');
end