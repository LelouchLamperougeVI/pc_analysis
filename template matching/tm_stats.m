function S=tm_stats(analysis_pre,analysis_post,sig)

if nargin<3
    sig=.99;
end

S=run_stats(analysis_pre,sig);
S=[S run_stats(analysis_post,sig)];

figure;
b=bar([S(1).fracEpochs; S(2).fracEpochs]');
b(1).FaceColor='r';
b(2).FaceColor='b';
title('fraction of significant epochs');
xlabel('compression factor');
ylabel('fraction epochs');
legend('pre-run','post-run');
xlim([0 length(S(1).fracEpochs)+1]);

figure;
boxscatter(S(1).sigCC,S(2).sigCC);


function boxscatter(data1,data2)
g=cellfun(@numel,data1);
g=arrayfun(@(x) -.1+x.*ones(1,g(x)),1:length(g),'uniformoutput',false);
groups=cell2mat(g);
data=cell2mat(data1');

g=cellfun(@numel,data2);
g=arrayfun(@(x) .1+x.*ones(1,g(x)),1:length(g),'uniformoutput',false);
groups=[groups cell2mat(g)];
data=[data; cell2mat(data2')];

groups=[groups (1:length(data1))];
data=[data; NaN(length(data1),1)];

groups=[groups (1:length(data1))+.5];
groups=[groups (1:length(data1))+.6];
data=[data; NaN(2*(length(data1)),1)];
boxplot(data,groups,'plotstyle','compact','color','rkbkk');
set(gca,'xtick',2:5:length(unique(groups)));
labels=num2cell(1:length(data1));
labels=cellfun(@num2str,labels,'uniformoutput',false);
set(gca,'xticklabel',labels);
title('correlation scores at significant epochs');
xlabel('compression factor');
ylabel('Spearman''s rho');
idx=findobj(gca,'Tag','Box');
legend(idx([5 3]),'pre-run','post-run');


function S=run_stats(analysis,sig)
C=analysis.C;
pval=analysis.C_pval;

idx=pval>sig;
S.fracEpochs=sum(idx)./size(pval,1); %fraction of significant epochs

S.sigCC=arrayfun(@(x) C(idx(:,x),x), 1:size(C,2), 'uniformoutput',false); %CCs during the significant epochs