function H=get_entropy(deconv,precision)
% helper function
% calculate entropy in each frame at specified precision

if islogical(deconv)
    p=sum(deconv,2);
    p=p./size(deconv,2);
    p=[p 1-p]';
    H=-sum(p.*log2(p),'omitnan');
else
    deconv=ca_filt(deconv);
    deconv=log(deconv);
    deconv(isinf(deconv))=nan;
    deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan');
    
    edges=[-inf min(deconv(:)):range(deconv(:))/precision:max(deconv(:))];
    edges(end)=inf;
    deconv(isnan(deconv))=min(min(deconv))-1;
    
    
    p=zeros(precision+1,size(deconv,1));
    for i=1:precision+1
        temp=deconv>=edges(i) & deconv<edges(i+1);
        p(i,:)=sum(temp,2)./size(deconv,1);
    end
    H=-sum(p.*log2(p),'omitnan');
end