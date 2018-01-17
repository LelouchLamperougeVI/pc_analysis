function [scores,pval,series]=pca_nw_msa(rest1,run,rest2,behavior,analysis)
% PCA-NW align

% [behavior,run]=convert_behavior(behavior,tcs,run);
% analysis=pc_batch_analysis(behavior,run);

% run=run(:,analysis.pc_list);
rest1=rest1(:,analysis.pc_list);
rest2=rest2(:,analysis.pc_list);

% runScore=pca_score(run,mean(diff(behavior.trials)));
[pc1,latent1]=pca_score(rest1,mean(diff(behavior.trials))/7);
[pc2,latent2]=pca_score(rest2,mean(diff(behavior.trials))/7);

pc1=pc1(:,latent1>1); % eigenvalue threshold set at 1 because zscored
pc2=pc2(:,latent2>1);

s1=pca_sequences(pc1,rest1);
% s2=get_run_sequence(behavior,tcs,run);
s3=pca_sequences(pc2,rest2);

s2=analysis.raw_stack(:,analysis.pc_list);
[~,s2]=max(s2);
[~,s2]=sort(s2);

%backwards replay
% s2=s2(end:-1:1);

ss2=mat2cell(s2,ones(1,size(s2,1)));

seq=horzcat(s1,ss2',s3);

[scores,pval]=mult_needleWunsch(seq);

pval=pval(length(s1)+1,:);
pval=find(pval<0.05);
series=[];
for i=1:length(pval)
    if pval(i)<=length(s1)
        series=[series get_template(s1{pval(i)},s2)];
    elseif pval(i)==length(s1)+1
        series=[series get_template(s2,s2)];
    else
        series=[series get_template(s3{pval(i)-length(s1)-1},s2)];
    end
end

function template=get_template(seq,match)

template=2.*double(match'==seq);
% template=2.*double(template==floor(template));

template=[template zeros(length(match),3) ones(length(match),2) zeros(length(match),3)];