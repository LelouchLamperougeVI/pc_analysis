function SI=get_si(raw_count,edges,Pi)

lambda=zeros(length(edges),length(raw_count),size(raw_count{1},2)); %P x bins x cells
for i=1:length(raw_count)
    temp=discretize(raw_count{i},edges);
    temp=temp+1;
    temp(isnan(temp))=1;
    temp=temp(:)+repelem((0:size(temp,2)-1)'.*length(edges),size(temp,1));
    temp=accumarray(temp,1,[length(edges)*size(raw_count{i},2) 1]);
    temp=reshape(temp,length(edges),size(raw_count{i},2));
%     temp=arrayfun(@(x) accumarray(temp(:,x),1,[length(edges) 1]),1:size(temp,2),'uniformoutput',false);
%     temp=cell2mat(temp);
%     temp=arrayfun(@(x) histcounts(raw_count{i}(:,x),edges,'binlimits',limits(:,x)), 1:size(raw_count{1},2),'uniformoutput',false);
%     temp=cell2mat(temp')';
    lambda(:,i,:)=permute(temp,[1 3 2]);
%     temp=sum(isnan(raw_count{i}) | isinf(raw_count{i}));
%     lambda(1,i,:)=temp;
end
lambda=lambda./sum(sum(lambda,1),2);
Pi=Pi./sum(Pi);
lambda_m=sum(lambda,2);
% Pi=sum(lambda,1);

SI=lambda.*log2(lambda./lambda_m./Pi);
SI(isnan(SI) | isinf(SI))=0;

SI=sum(sum(SI,2),1);
SI=reshape(SI,1,[]);

H_x=lambda_m.*log2(lambda_m);
H_x(isinf(H_x) | isnan(H_x))=0;
H_x=-sum(H_x);
H_x=shiftdim(H_x,1);
SI=SI./H_x;