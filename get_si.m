function SI=get_si(raw_count,edges,Pi)

% lambda=sum(raw_psth>0)./size(raw_psth,1); %%%%%%=>lambda_test=sum(raw_psth>0)./sum(sum(raw_psth>0)) I think that this is a conditional propability!!!!!
% SI=Pi.*lambda./(sum(lambda,2)./size(raw_psth,2)).*log2(lambda./(sum(lambda,2)./size(raw_psth,2)));
% SI(isnan(SI))=0;
% SI=sum(SI,2);
% SI=reshape(SI,1,[]);

% MI version
% maxi=cellfun(@max, raw_count,'uniformoutput',false);
% maxi=cell2mat(maxi');
% maxi=max(maxi);
% mini=cellfun(@min, raw_count,'uniformoutput',false);
% mini=cell2mat(mini');
% mini=min(mini);
% limits=[mini;maxi];
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

% SI(1,:,:)=[];

SI=sum(sum(SI,2),1);
SI=reshape(SI,1,[]);
% temp=sum(raw_count);
% P1=temp./Pi;
% P0=1-P1;
% Pi=Pi./sum(Pi);
% temp=[Pi.*P1;Pi.*P0];
% % temp=Pi.*P1;
% temp=temp.*log2(temp./sum(temp,1)./sum(temp,2));
% temp(isnan(temp))=0;
% SI=sum(sum(temp));
% % SI=sum(temp);
% SI=reshape(SI,1,[]);

% % % paper version
% lambda_m=sum(sum(raw_count))./sum(Pi);
% lambda=sum(raw_count)./Pi;
% Pi=Pi./sum(Pi);
% SI=lambda./lambda_m.*Pi.*log2(lambda./lambda_m);
% % SI=lambda.*Pi.*log2(lambda./lambda_m);
% SI(isnan(SI))=0;
% SI=sum(SI);
% SI=reshape(SI,1,[]);