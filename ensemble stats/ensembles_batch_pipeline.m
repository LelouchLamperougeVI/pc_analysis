%% Batch analyze ensembles rest-run-rest

precision=5;

d=dir('2017*');

for i=1:length(d)
    cd(d(i).name);
    
    sessions=dir;
    sessions=arrayfun(@(x) str2double(sessions(x).name), 1:length(sessions));
    sessions=mod(sessions,3);
    sessions(isnan(sessions))=1;
    sessions=find(~sessions)-3; %list of behavior folders
    
    for j=1:length(sessions)
        try
        load([num2str(sessions(j)) '\behavior.mat']);
        load([num2str(sessions(j)) '\Plane1\deconv.mat']);
        load([num2str(sessions(j)) '\Plane1\timecourses.mat']);
        
        [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
        analysis=pc_batch_analysis(behavior,deconv,'test','mixed','bins',80);
        order=get_order(analysis);
        
        assemblies=cluster_mi(deconv,'prune',10,'plotFlag',false,'shuffle',20,'sig',5,'precision',precision);
        sce=knn_sce(deconv,assemblies);
        
        save([num2str(sessions(j)) '\react.mat'],'order','analysis','deconv','assemblies','sce');
        
        %pre-rest
        load([num2str(sessions(j)-1) '\behavior.mat']);
        load([num2str(sessions(j)-1) '\Plane1\deconv.mat']);
        
        vel=behavior.speed_raw;
        thres=noRun(vel);
        deconv(vel>thres,:)=[];
        
        assemblies=cluster_mi(deconv,'prune',10,'plotFlag',false,'shuffle',20,'sig',5,'precision',precision);
        sce=knn_sce(deconv,assemblies);
        
        save([num2str(sessions(j)-1) '\react.mat'],'order','deconv','assemblies','sce');
        
        %post_rest
        load([num2str(sessions(j)+1) '\behavior.mat']);
        load([num2str(sessions(j)+1) '\Plane1\deconv.mat']);
        
        vel=behavior.speed_raw;
        thres=noRun(vel);
        deconv(vel>thres,:)=[];
        
        assemblies=cluster_mi(deconv,'prune',10,'plotFlag',false,'shuffle',20,'sig',5,'precision',precision);
        sce=knn_sce(deconv,assemblies);
        
        save([num2str(sessions(j)+1) '\react.mat'],'order','deconv','assemblies','sce');
        catch
        end
    end
    
    cd ..
end

%% Stats

windows=[]; %window number
prc_ass_ovl=[]; %fraction overlap with run assemblies
prc_pc_ovl=[]; %fraction overlap with place cells
prc_ass_null=[]; %fraction overlap with run by chance
prc_pc_null=[]; %fraction overlap with place cells by chance
num_ens=[]; %number of ensembles
ens_size=[]; %size of ensembles
num_events=[]; %number of sce events - reported as percentage of total recording
pc_frac=[];
days=[];

% fid=fopen('rsc_log.csv','r');
fid=fopen('rsc37_log.csv','r');
tb=textscan(fid,'%s');
fclose(fid);
date=textscan(tb{1}{1},'%{dd/MM/uuuu}D','delimiter',','); date=date{1};
win=textscan(tb{1}{2},'%*3q%1d"','delimiter',','); win=win{1};
sess=textscan(tb{1}{3},'"%1d"','delimiter',','); sess=sess{1};

while length(win)<length(date); win=[win;0]; end


d=dir('2017*');

for i=1:length(d)
    cd(d(i).name);
    
    sessions=dir;
    sessions=arrayfun(@(x) str2double(sessions(x).name), 1:length(sessions));
    sessions=mod(sessions,3);
    sessions(isnan(sessions))=1;
    sessions=find(~sessions)-3; %list of behavior folders
    
    for j=1:length(sessions)
        if exist([num2str(sessions(j)) '\react.mat'],'file') && exist([num2str(sessions(j)-1) '\react.mat'],'file') && exist([num2str(sessions(j)+1) '\react.mat'],'file')
            idx=datetime(d(i).name,'inputformat','uuuu_MM_dd')==date & ceil(sessions(j)/3)==sess;
            days=[days;date(idx)];
            windows=[windows;win(idx)];
            
            load([num2str(sessions(j)-1) '\react.mat']);
            assemblies1=assemblies;
            sce1=sce;
            deconv1=deconv;
            
            load([num2str(sessions(j)) '\react.mat']);
            assemblies2=assemblies; 
            sce2=sce;
            deconv2=deconv;
            
            load([num2str(sessions(j)+1) '\react.mat']);
            assemblies3=assemblies;
            sce3=sce;
            deconv3=deconv;
            
            pc_list=analysis.pc_list;
            
            temp1=0;temp2=0;
            if ~isempty(sce1)
                for k=1:length(sce1)
                    temp=(deconv1(:,assemblies1{k})>0) & sce1{k};
                    temp=logical(sum(temp,2));
                    temp1=temp+temp1;
                end
                temp1=sum(logical(temp1))/length(temp1);
            end
            if ~isempty(sce3)
                for k=1:length(sce3)
                    temp=(deconv3(:,assemblies3{k})>0) & sce3{k};
                    temp=logical(sum(temp,2));
                    temp2=temp+temp2;
                end
                temp2=sum(logical(temp2))/length(temp2);
            end
            num_events=[num_events; temp1 temp2];
            
            num_ens=[num_ens; length(assemblies1) length(assemblies3)];
            assemblies1=cell2mat(assemblies1);
            assemblies2=cell2mat(assemblies2);
            assemblies3=cell2mat(assemblies3);
            ens_size=[ens_size; length(assemblies1)./size(deconv,2) length(assemblies3)./size(deconv,2)];
            pc_frac=[pc_frac; length(pc_list)./size(deconv,2)];
            
            prc_ass_ovl=[prc_ass_ovl; (length(assemblies1) + length(assemblies2) - length(setxor(assemblies1,assemblies2))) / 2 / (length(assemblies1) + length(assemblies2))...
                                       (length(assemblies3) + length(assemblies2) - length(setxor(assemblies3,assemblies2))) / 2 / (length(assemblies3) + length(assemblies2)) ];
            prc_ass_null=[prc_ass_null; (length(assemblies2)/size(deconv,2)) * (length(assemblies1)/size(deconv,2)) ...
                                        (length(assemblies2)/size(deconv,2)) * (length(assemblies3)/size(deconv,2))];
                                    
            prc_pc_ovl=[prc_pc_ovl; (length(assemblies1) + length(pc_list) - length(setxor(assemblies1,pc_list))) / 2 / length(assemblies1)...
                                       (length(assemblies3) + length(pc_list) - length(setxor(assemblies3,pc_list))) / 2 / length(assemblies3) ];
            prc_pc_null=[prc_pc_null; (length(pc_list)/size(deconv,2)) ...
                                        (length(pc_list)/size(deconv,2))];
            
            
        end
    end
    
    cd ..;
end

%%

g={'RSC','M1/2','S1', 'PPC'};
windows(windows==1 | windows==2)=1;
windows(windows==3 | windows==4 | windows==8)=2; 
windows(windows==6 | windows==7)=3; 
windows(windows==5)=4;

%% PC overlap
figure;
subplot(1,2,1);hold on;
boxplot(prc_pc_ovl(:,1),windows,'colors','r');
boxplot(prc_pc_null(:,1),windows,'colors','k');
title('Pre-rest');
ylim([0 1]);
ylabel('fraction');
xticklabels(g);

subplot(1,2,2);hold on;
boxplot(prc_pc_ovl(:,2),windows,'colors','r');
boxplot(prc_pc_null(:,2),windows,'colors','k');
title('Post-rest');
ylim([0 1]);
ylabel('fraction');
xticklabels(g);

plot(nan,'r');
plot(nan,'k');
legend('actual','by chance');

figure;
subplot(1,2,1);hold on;
boxplot(prc_pc_ovl(:,1)-prc_pc_null(:,1),windows);
title('Pre-rest');
ylabel('\Delta fraction');
xticklabels(g);

subplot(1,2,2);hold on;
boxplot(prc_pc_ovl(:,2)-prc_pc_null(:,2),windows);
title('Post-rest');
ylabel('\Delta fraction');
xticklabels(g);


%% Ensembles overlap
figure;
subplot(1,2,1);hold on;
boxplot(prc_ass_ovl(:,1),windows,'colors','r');
boxplot(prc_ass_null(:,1),windows,'colors','k');
title('Pre-rest');
ylim([0 1]);
ylabel('fraction');
xticklabels(g);

subplot(1,2,2);hold on;
boxplot(prc_ass_ovl(:,2),windows,'colors','r');
boxplot(prc_ass_null(:,2),windows,'colors','k');
title('Post-rest');
ylim([0 1]);
ylabel('fraction');
xticklabels(g);

plot(nan,'r');
plot(nan,'k');
legend('actual','by chance');

%% Ensembles sizes
figure; hold on
boxplot(ens_size(:,1),windows,'colors','r');
boxplot(ens_size(:,2),windows,'colors','k');
ylabel('fraction cells');
xticklabels(g);

plot(nan,'r');
plot(nan,'k');
legend('pre','post');


%% percentage events
figure;
subplot(1,2,1);
boxplot(num_events(:,1).*100,windows,'colors','r');
ylabel('ensembles mean firing rate (%)');
title('Pre-rest');
xticklabels(g);
subplot(1,2,2);
boxplot(num_events(:,2).*100,windows,'colors','k');
title('Pre-rest');
xticklabels(g);

% plot(nan,'r');
% plot(nan,'k');
% legend('pre','post');

%%
figure;
bins=0:max(num_ens(:));

pre=zeros(length(bins),length(g));
post=pre;

for i=1:length(g)
    idx=num_ens(windows==i,:)==permute(bins,[3 1 2]);
    idx=squeeze(sum(idx))';
    
    pre(:,i)=idx(:,1);
    post(:,i)=idx(:,2);
end

subplot(1,2,1);hold on;
bar(pre);
title('Pre-rest');
ylabel('# of occurrences');
xlabel('# of ensembles');

subplot(1,2,2);hold on;
bar(post);
title('Post-rest');
ylabel('# of occurrences');
xlabel('# of ensembles');

legend(g);