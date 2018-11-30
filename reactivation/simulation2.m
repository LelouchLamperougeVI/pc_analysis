m=12000;n=200;
ass=zeros(3,n);
ass(1,1:50)=1;
ass(2,41:60)=1;
ass(3,58:62)=1;

jitter=4;

q=bin_poisson_sim(m,n,'R',.05,'assemblies',ass,'prob',.05,'miss',[10 4 1],'jitter',jitter);

figure;
subplot(1,2,1);
imagesc(corr(q));
axis square
title('Pearson Corr');

subplot(1,2,2);
imagesc(get_mi(q));
axis square
title('Mutual Info');

%% Accuracy as function of jitter (normal)
m=12000;n=200;
ass=zeros(3,n);
ass(1,1:50)=1;
ass(2,41:60)=1;
ass(3,58:62)=1;

true_ass=1:62;

precision=5;

it=10; %iterations per sd
sd=0:10;

acc=zeros(it,length(sd),3); % 5%, 1%, lopes
fp=acc;
tp=acc;
% detected_size=acc;
% e_size=acc;

for j=1:length(sd)
    for i=1:it
        q=bin_poisson_sim(m,n,'R',.05,'assemblies',ass,'prob',.05,'miss',[10 4 1],'jitter',sd(j),'jitter_type','uniform');
        q=logical(q);
        
        assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',100,'sig',5,'precision',precision);
        assemblies=unique(cell2mat(assemblies));
        tp(i,j,1)=sum(ismember(assemblies,true_ass));
        fp(i,j,1)=sum(ismember(assemblies,setxor(true_ass,1:n)));
        acc(i,j,1)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                    (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        
        assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',100,'sig',1,'precision',precision);
        assemblies=unique(cell2mat(assemblies));
        tp(i,j,2)=sum(ismember(assemblies,true_ass));
        fp(i,j,2)=sum(ismember(assemblies,setxor(true_ass,1:n)));
        acc(i,j,2)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                    (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        
        try
            assemblies=lopes_pca(q,0,false);
            assemblies=unique(cell2mat(assemblies));
            tp(i,j,3)=sum(ismember(assemblies,true_ass));
            fp(i,j,3)=sum(ismember(assemblies,setxor(true_ass,1:n)));
            acc(i,j,3)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                        (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        catch
        end
    end
    j/length(sd)
end

figure;
errorbar(sd,mean(acc(:,:,1)),sem(acc(:,:,1)),'-b')
hold on
errorbar(sd,mean(acc(:,:,2)),sem(acc(:,:,2)),'-r')
errorbar(sd,mean(acc(:,:,3)),sem(acc(:,:,3)),'-k')
legend('MI 5%','MI 1%','Lopes PCA')
xlabel('jitter (\pm samples; uniform noise)')
ylabel('accuracy (fraction)')

figure;
errorbar(sd,mean(tp(:,:,1)),sem(tp(:,:,1)),'-b')
hold on
errorbar(sd,mean(tp(:,:,2)),sem(tp(:,:,2)),'-r')
errorbar(sd,mean(tp(:,:,3)),sem(tp(:,:,3)),'-k')
errorbar(sd,mean(fp(:,:,1)),sem(fp(:,:,1)),'--b')
errorbar(sd,mean(fp(:,:,2)),sem(fp(:,:,2)),'--r')
errorbar(sd,mean(tp(:,:,3)),sem(tp(:,:,3)),'--k')
legend('MI 5%','MI 1%','Lopes PCA')
xlabel('jitter (\pm samples; uniform noise)')
ylabel('true positives (count; solid lines)')
yl=ylim;
yyaxis right
ylabel('false positives (count; dotted lines)')
ylim(yl);


%% Background fr vs ensembles fr
m=12000;n=200;
ass=zeros(3,n);
ass(1,1:50)=1;
ass(2,41:60)=1;
ass(3,58:62)=1;

true_ass=1:62;

precision=5;

bg_fr=.01:.01:.1;
e_fr=.01:.01:.1;
jitter=0:2;

acc=zeros(length(bg_fr),length(e_fr),length(jitter),3); % 5%, 1%, lopes

for j=1:length(bg_fr)
    for i=1:length(e_fr)
        for k=1:length(jitter)
            q=bin_poisson_sim(m,n,'R',bg_fr(j),'assemblies',ass,'prob',e_fr(i),'miss',[10 4 1],'jitter',jitter(k),'jitter_type','uniform');
            q=logical(q);

            assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',100,'sig',5,'precision',precision);
            assemblies=unique(cell2mat(assemblies));
            acc(j,i,k,1)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                        (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));

            assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',100,'sig',1,'precision',precision);
            assemblies=unique(cell2mat(assemblies));
            acc(j,i,k,2)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                        (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));

            try
                assemblies=lopes_pca(q,0,false);
                assemblies=unique(cell2mat(assemblies));
                acc(j,i,k,3)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                            (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
            catch
            end
        end
    end
    j/length(bg_fr)
end


figure;
for i=1:3
    for j=1:3
        subplot(3,3,3*(i-1)+j);
        imagesc(e_fr,bg_fr,acc(:,:,i,j));
        if i==1
            switch j
                case 1
                    title('MI 5%');
                case 2
                    title('MI 1%');
                case 3
                    title('Lopes PCA');
            end
        end
        if i==3 && j==2
            xlabel('ensembles firing rate');
        end
        if i==2 && j==1
            ylabel('background firing rate');
        end
        axis square;
        c=colorbar;
        caxis([0 1]);
        if i==2 && j==3
            c.Label.String='accuracy';
        end
    end
end



%%
m=12000;n=200;
precision=5;
it=10; %iterations per sd
s=2:10;
acc=zeros(it,length(s),2);
detected_size=acc;

for j=s
    e_size=j^2;
    ass=zeros(1,n);
    ass(1:e_size)=1;
    
    true_ass=1:e_size;
    for i=1:it
        q=bin_poisson_sim(m,n,'R',.05,'assemblies',ass,'prob',.05,'miss',round(.2*e_size),'jitter',1);
        assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',10,'sig',5,'precision',precision);
        assemblies=unique(cell2mat(assemblies));
        detected_size(i,j-1,2)=length(assemblies)/length(true_ass);
        acc(i,j-1,2)=1 - length(setxor(true_ass,assemblies))/(length(true_ass)+length(assemblies));
        assemblies=lopes_pca(q,1,false);
        assemblies=unique(cell2mat(assemblies));
        detected_size(i,j-1,1)=length(assemblies)/length(true_ass);
        acc(i,j-1,1)=1 - length(setxor(true_ass,assemblies))/(length(true_ass)+length(assemblies));
    end
end


%% Separability defined as t-statistic as a function of jitter and ensemble size
m=12000;n=200;
s=1:10:n-1;
sd=0:10;
tstat=zeros(length(s),length(sd),2);
pval=tstat;

for i=1:length(s)
    e_size=s(i);
    ass=zeros(1,n);
    ass(1:e_size)=1;
    
    for j=1:length(sd)
        q=bin_poisson_sim(m,n,'R',.05,'assemblies',ass,'prob',.05,'miss',round(.2*e_size),'jitter',sd(j));
        
        R=1-corr(q);
        I=get_mi(q);
        
        x=[squareform(R(1:e_size,1:e_size));squareform(I(1:e_size,1:e_size))];
        y=[squareform(R(e_size+1:end,e_size+1:end));squareform(I(e_size+1:end,e_size+1:end))];
        
        [~,pval(i,j,1),~,stats]=ttest2(x(1,:),y(1,:));
        tstat(i,j,1)=stats.tstat;
        [~,pval(i,j,2),~,stats]=ttest2(x(2,:),y(2,:));
        tstat(i,j,2)=stats.tstat;
    end
    s(i)
end


%% Silhouette scores for place cells classification
list=dir('*tread*new.mat');

s_mi=[];s_corr=[];
for i=1:length(list)
    load(list(i).name);
    deconv=analysis.deconv;
    clust=zeros(1,size(deconv,2));
    clust(analysis.pc_list)=1;
    [~,D_mi]=knn_mi(deconv);
    Y_mi=cmdscale(D_mi);
    D_corr=corr(deconv);
    D_corr=1-abs(D_corr);
    Y_corr=cmdscale(D_corr);
    
    s=silhouette(Y_mi,clust);
    s_mi=[s_mi;s(analysis.pc_list)];
    s=silhouette(Y_corr,clust);
    s_corr=[s_corr;s(analysis.pc_list)];
end

figure
[f,x]=ecdf(s_mi);
plot(x,f,'r');
hold on
[f,x]=ecdf(s_corr);
plot(x,f,'k');
xlabel('silhouette score');
ylabel('cummulative frequency');
plot([0 0],[0 1],'--b');
legend('MI','Pearson');
[~,p]=kstest2(s_mi,s_corr,'tail','smaller');
title(['one-tailed KS test p = ' num2str(p)]);

%% MI vs Pearson; correlation between spikes & jitter

m=12000;n=200;
ass=zeros(3,n);
ass(1,1:50)=1;
ass(2,41:60)=1;
ass(3,58:62)=1;

true_ass=1:62;
false_ass=setxor(true_ass,1:n);

precision=5;

jitter=0:10;
r=.05:.05:.95;

acc=zeros(length(jitter),length(r),3); %mi,pearson,lopes
for i=1:length(jitter)
    for j=1:length(r)
        q=bin_poisson_sim(m,n,'R',.05,'assemblies',ass,'prob',.05,'miss',[10 4 1],'jitter',jitter(i),'corr',r(j),'jitter_type','uniform','seed',deconv);
        
        assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',10,'sig',5,'precision',precision,'method','mi','nofilt',false);
        assemblies=unique(cell2mat(assemblies));
        acc(i,j,1)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                    (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        
        assemblies=cluster_mi(q,'prune',4,'plotFlag',false,'shuffle',10,'sig',5,'precision',precision,'method','pearson','nofilt',false);
        assemblies=unique(cell2mat(assemblies));
        acc(i,j,2)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                    (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        
        try
            assemblies=lopes_pca(q,0,false);
            assemblies=unique(cell2mat(assemblies));
            acc(i,j,3)=(sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass))) / ...
                        (sum(ismember(assemblies,true_ass)) + sum(ismember(setxor(assemblies,1:n),false_ass)) + sum(ismember(assemblies,false_ass)) + sum(ismember(setxor(assemblies,1:n),true_ass)));
        catch
        end
    end
end


figure;
subplot(1,3,1);
imagesc(r,jitter,acc(:,:,1));
axis square
caxis([0 1]);
title('MI');
ylabel('jitter (\pm samples; uniform distribution)');
colorbar;

subplot(1,3,2);
imagesc(r,jitter,acc(:,:,2));
axis square
caxis([0 1]);
title('Pearson');
xlabel('spikes correlation');
colorbar;

subplot(1,3,3);
imagesc(r,jitter,acc(:,:,3));
axis square
caxis([0 1]);
title('Lopes');
c=colorbar;
c.Label.String='accuracy';