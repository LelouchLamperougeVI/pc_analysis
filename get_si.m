function SI=get_si(raw_count,raw_psth,Pi)

% lambda=sum(raw_psth>0)./size(raw_psth,1); %%%%%%=>lambda_test=sum(raw_psth>0)./sum(sum(raw_psth>0)) I think that this is a conditional propability!!!!!
% SI=Pi.*lambda./(sum(lambda,2)./size(raw_psth,2)).*log2(lambda./(sum(lambda,2)./size(raw_psth,2)));
% SI(isnan(SI))=0;
% SI=sum(SI,2);
% SI=reshape(SI,1,[]);

% MI version
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

% % paper version
lambda_m=sum(sum(raw_count))./sum(Pi);
lambda=sum(raw_count)./Pi;
Pi=Pi./sum(Pi);
SI=lambda./lambda_m.*Pi.*log2(lambda./lambda_m);
% SI=lambda.*Pi.*log2(lambda./lambda_m);
SI(isnan(SI))=0;
SI=sum(SI);
SI=reshape(SI,1,[]);