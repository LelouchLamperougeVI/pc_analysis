function te=embed_dimensions_gpu(x,y,s,l,order,chunked)
% x -> y embedding

idx=[y+s.*(0:order) x+s.*(1:order)];
p_joint=accumarray(chunked(:,idx),1)./l;

idx=[y+s.*(1:order) x+s.*(1:order)];
p_cond_xy=accumarray(chunked(:,idx),1)./l;
p_cond_xy=permute(p_cond_xy,[ndims(p_cond_xy)+1 1:ndims(p_cond_xy)]); %inject singletons
p_cond_xy=p_joint./p_cond_xy;

idx=y+s.*(0:order);
p_joint_yy=accumarray(chunked(:,idx),1)./l;

idx=y+s.*(1:order);
p_joint_y=accumarray(chunked(:,idx),1)./l;
p_joint_y=permute(p_joint_y,[ndims(p_joint_y)+1 1:ndims(p_joint_y)]);
p_cond_y=p_joint_yy./p_joint_y;

te=p_joint.*log2(p_cond_xy./p_cond_y);
te=sum(te(:),'omitnan');
end