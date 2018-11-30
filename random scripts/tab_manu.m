figure


subplot(1,3,1);
[f,x]=ecdf(tab_SI);
disp(['tab: ' num2str(mean(tab_SI)) ' +- ' num2str(sem(tab_SI))]);
plot(x,f,'r');
[f,x]=ecdf(tread_SI);
disp(['tread: ' num2str(mean(tread_SI)) ' +- ' num2str(sem(tread_SI))]);
hold on
plot(x,f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('mutual info. (bits)');
ylabel('cum. prob.');
[~,p]=ttest2(tab_SI,tread_SI)

subplot(1,3,2);
[f,x]=ecdf(tab_sparsity);
disp(['tab: ' num2str(mean(tab_sparsity)) ' +- ' num2str(sem(tab_sparsity))]);
plot(x,f,'r');
[f,x]=ecdf(tread_sparsity);
disp(['tread: ' num2str(mean(tread_sparsity)) ' +- ' num2str(sem(tread_sparsity))]);
hold on
plot(x,f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('sparsity (%)');
ylabel('cum. prob.');
[~,p]=ttest2(tab_sparsity,tread_sparsity)

subplot(1,3,3);
idx1=cellfun(@isempty, tab_width);
idx1=cell2mat(tab_width(~idx1)');
[f,x]=ecdf(idx1(:,1));
disp(['tab: ' num2str(mean(idx1(:,1))) ' +- ' num2str(sem(idx1(:,1)))]);
plot(x.*(100/80),f,'r');
idx2=cellfun(@isempty, tread_width);
idx2=cell2mat(tread_width(~idx2)');
[f,x]=ecdf(idx2(:,1));
disp(['tread: ' num2str(mean(idx2(:,1))) ' +- ' num2str(sem(idx2(:,1)))]);
hold on
plot(x.*(100/80),f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('place field width (cm)');
ylabel('cum. prob.');
[~,p]=ttest2(idx1(:,1),idx2(:,1))


%% behavior start, middle, end
behavior=analysis.behavior;
bins=length(analysis.Pi);

vel_thres=noRun(behavior.unit_vel);
vel_thres=behavior.unit_vel>vel_thres | behavior.unit_vel<-vel_thres;
behavior.unit_vel=behavior.unit_vel(vel_thres);
behavior.unit_pos=behavior.unit_pos(vel_thres);
behavior.frame_ts=behavior.frame_ts(vel_thres);

vel=zeros(bins,length(behavior.trials)-1);
idx=linspace(min(behavior.unit_pos),max(behavior.unit_pos),bins+1);
idx=linspace(-100,0,bins+1);
for i=1:length(behavior.trials)-1
    for j=1:length(idx)-1
        tmp=behavior.unit_pos>idx(j) & behavior.unit_pos<=idx(j+1) & behavior.frame_ts>behavior.trials(i) & behavior.frame_ts<=behavior.trials(i+1);
        vel(j,i)=mean(behavior.unit_vel(tmp));
    end
end
sd=4/analysis.vr_length*bins;

try
    idx=get_head(flipud(isnan(vel)));
    idx=flipud(idx);
    idx=arrayfun(@(x) find(idx(:,x),1), 1:size(idx,2));
    for i=1:length(idx)
        vel(1:idx(i),i)=0;
    end
catch
end
vel=fillmissing(vel,'linear');
vel=fast_smooth(vel,sd);

start=mean(vel(1:bins*.05,:));
middle=mean(vel(.25*bins:.75*bins,:));
ending=mean(vel(end-bins*.05+1:end,:));

errorbar(mean([start' middle' ending']),sem([start' middle' ending']));
hold on

%%
xticklabels({'start','','middle','','end'});
ylabel('mean velocity (cm/s)');
xlim([.5 3.5]);
ylim([0 30]);
legend(split(num2str(1:6),'  '))


%% 2d navigation stats

list=dir;

rp_ratio=zeros(length(list)-2,1);
day=rp_ratio;
numRew=rp_ratio;
numPun=rp_ratio;
numTrials=rp_ratio;
for i=3:length(list)
    behavior=import_2d(list(i).name,list(i).folder);
    save([list(i).name(1:end-4) '_beh'],'behavior');
    
    numTrials(i-2)=length(behavior.trial)-1;
    
    reward=behavior.reward(behavior.reward>behavior.trial(1) & behavior.reward<=behavior.trial(end));
    punish=behavior.punish(behavior.punish>behavior.trial(1) & behavior.punish<=behavior.trial(end));
    
    numRew(i-2)=length(reward);
    numPun(i-2)=length(punish);
    
    rp_ratio(i-2)=numRew(i-2)/numPun(i-2);
    
    day(i-2)=mod(behavior.date,100);
end

save('analysis','day','numPun','numRew','numTrials','rp_ratio');

%% Rerun stats for last 10 min and last 20 trials
list=dir('P*');
day=cell(1,4);
numRew=day;
numPun=day;
numTrials=day;
x=day;
y=day;
t=day;
ts=day;

for i=1:length(list)
    idx=dir([list(i).name '\*_beh.mat']);
    day{i}=zeros(length(idx),1);
    numRew{i}=day{i};
    numPun{i}=day{i};
    numTrials{i}=day{i};
    x{i}=cell(1,length(idx));
    y{i}=cell(1,length(idx));
    ts{i}=cell(1,length(idx));
    t{i}=day{i};
    for j=1:length(idx)
        load([list(i).name '\' idx(j).name]);
        day{i}(j)=mod(behavior.date,100);
        
        tf=[behavior.ts(behavior.trial(1)) behavior.ts(behavior.trial(end))]; %time-frame for last 10 min of recording
        tf=(diff(tf)/60-10)*60 + tf(1);
        tf=behavior.trial(find(behavior.ts(behavior.trial)>tf,1))-1;
        
        numRew{i}(j)=sum(behavior.reward>tf & behavior.reward<=behavior.trial(end));
        numPun{i}(j)=sum(behavior.punish>tf & behavior.punish<=behavior.trial(end));
        numTrials{i}(j)=sum(behavior.trial>tf);
        
        x{i}{j}=behavior.x(behavior.trial(end-19):behavior.trial(end));
        y{i}{j}=behavior.y(behavior.trial(end-19):behavior.trial(end));
        ts{i}{j}=behavior.ts(behavior.trial(end-19):behavior.trial(end));
        
        t{i}(j)=(behavior.ts(behavior.trial(end))-behavior.ts(tf))/60;
    end
end

%%
rew=[numRew{1}([1 4])./numTrials{1}([1 4]) numRew{2}([1 4])./numTrials{2}([1 4]) numRew{3}([1 4])./numTrials{3}([1 4]) numRew{4}([1 4])./numTrials{4}([1 4])];
pun=[numPun{1}([1 4])./numTrials{1}([1 4]) numPun{2}([1 4])./numTrials{2}([1 4]) numPun{3}([1 4])./numTrials{3}([1 4]) numPun{4}([1 4])./numTrials{4}([1 4])];

figure;
plot([rew(:,[3 4 1 2]); NaN(1,4); pun(:,[3 4 1 2])]);
xlim([0 6]);
ylim([.5 2]);
xticklabels({'','day 1','day4','','day 1','day4'});


t_per_min=[numTrials{1}([1 4])./t{1}([1 4]) numTrials{2}([1 4])./t{2}([1 4]) numTrials{3}([1 4])./t{3}([1 4]) numTrials{4}([1 4])./t{4}([1 4])];
figure;
plot(t_per_min(:,[3 4 1 2]));
xlim([0 3]);
ylim([0 10]);
xticklabels({'','','day 1','','day4'});


xn=20; %10 cm bins
yn=14;
xedges=linspace(-100,100,xn);
yedges=linspace(0,140,yn);
c=zeros(2,4);
avg_vel=c;

for i=1:4
    for j=1:2
        xx=x{i}{j^2};
        yy=y{i}{j^2};
        
        vel=diff([xx yy]);
        vel=sqrt(sum(vel.^2,2))./diff(ts{i}{j^2});
        vel(isoutlier(vel))=0;
        thres=noRun(vel);
        avg_vel(j,i)=mean(vel(vel>thres));
        
        vel=[0; vel];
        
        xx(vel<thres)=[];
        yy(vel<thres)=[];
        
        p=histcounts2(xx,yy,xedges,yedges,'normalization','probability');
        temp=p(:).*log(p(:));
        temp(isnan(temp))=0;
        c(j,i)=-sum(temp)/log(xn*yn);
    end
end
figure;
plot(c(:,[3 4 1 2]));
xlim([0 3]);
ylim([0 1]);
xticklabels({'','','day 1','','day4'});
figure;
plot(avg_vel(:,[3 4 1 2]));
xlim([0 3]);
ylim([0 20]);
xticklabels({'','','day 1','','day4'});

%%
list=dir;

rew=[];
pun=[];
day=[];

count=1;
for i=4:2:length(list)
    load(list(i).name);
    for j=1:length(behavior.trial)-1
        idx=behavior.reward(behavior.reward>behavior.trial(j) & behavior.reward<=behavior.trial(j+1));
        rew=[rew length(idx)];
        idx=behavior.punish(behavior.punish>behavior.trial(j) & behavior.punish<=behavior.trial(j+1));
        pun=[pun length(idx)];
        day=[day count];
    end
    count=count+1;
end


%%
bar(pool([1 4],[3 4 1 2])')
legend("left correct", "left incorrect")
ylim([0 .25])

%%
figure
bar(occ_pool([1 2],:)');
legend("left rew trials", "right rew trials");
ylim([0 .3])
figure
bar(occ_pool([3 4],:)');
legend("right rew trials", "left rew trials");
ylim([0 .3])
figure
bar(rew_pool([1 2],:)');
legend("left rewards", "left punishments");
ylim([0 1.5])
figure
bar(rew_pool([3 4],:)');
legend("right rewards", "right punishments");
ylim([0 1.5])


%% Trajectory density/Concentration
xn=20; %10 cm bins
yn=14;
xedges=linspace(-100,100,xn);
yedges=linspace(0,140,yn);

p=histcounts2(behavior.x(behavior.trial(1):behavior.trial(end)),behavior.y(behavior.trial(1):behavior.trial(end)),xedges,yedges,'normalization','probability');
c=p(:).*log(p(:));
c(isnan(c))=0;
c=-sum(c)/log(xn*yn);

%% rew/pun per trial
figure; hold on;

plot([numRew(1)./numTrials(1) numRew(4)./numTrials(4) nan numPun(1)./numTrials(1) numPun(4)./numTrials(4)])

plot([numRew(1)./numTrials(1) numRew(2)./numTrials(2) numRew(3)./numTrials(3) numRew(4)./numTrials(4) nan...
    numPun(1)./numTrials(1) numPun(2)./numTrials(2) numPun(3)./numTrials(3) numPun(4)./numTrials(4)])

xlim([0 6]);
ylim([.5 2]);
xticklabels({'','day 1','day4','','day 1','day4'});

%% trials/min
figure; hold on;

plot([t_per_min(1) t_per_min(4)]);
xlim([0 3]);
ylim([0 10]);
xticklabels({'','','day 1','','day4'});

%% trials concentration
figure; hold on;

plot([cct(1) cct(4)]);
xlim([0 3]);
ylim([.5 1]);
xticklabels({'','','day 1','','day4'});


%%
list=dir('VR*.mat');
rew_per_trial=cell(1,4);
pun_per_trial=cell(1,4);
trials_duration=cell(1,4);
for i=1:4
    load(list(i).name);
    rew_per_trial{i}=arrayfun(@(x) sum( behavior.reward>behavior.trial(x) & behavior.reward<=behavior.trial(x+1) ), 1:length(behavior.trial)-1);
    pun_per_trial{i}=arrayfun(@(x) sum( behavior.punish>behavior.trial(x) & behavior.punish<=behavior.trial(x+1) ), 1:length(behavior.trial)-1);
    trials_duration{i}=arrayfun(@(x) behavior.ts(behavior.trial(x+1)) - behavior.ts(behavior.trial(x)), 1:length(behavior.trial)-1);
end

rew_per_trial=cell2mat(rew_per_trial);
pun_per_trial=cell2mat(pun_per_trial);
trials_duration=cell2mat(trials_duration);

bins=20;
rew_per_trial=arrayfun(@(x) mean(rew_per_trial(x:x+bins-1)), 1:bins:length(rew_per_trial)-bins);
pun_per_trial=arrayfun(@(x) mean(pun_per_trial(x:x+bins-1)), 1:bins:length(pun_per_trial)-bins);
trials_duration=arrayfun(@(x) mean(trials_duration(x:x+bins-1)), 1:bins:length(trials_duration)-bins);


%%
list=dir('P*');
figure; hold on;
for i=1:length(list)
    load([list(i).name '\analysis.mat']);
    numRew(5)=[];
    numTrials(5)=[];
    plot(numRew./numTrials);
end


%% Redo pc analysis
list=dir('*_tab_*_new.mat');

tab_SI=[];
tab_sparsity=[];
tab_width=[];
for i=1:length(list)
    load(list(i).name);
    
    tab_SI=[tab_SI analysis.SI];
    tab_sparsity=[tab_sparsity analysis.sparsity];
    tab_width=[tab_width analysis.width(analysis.pc_list)];
end

list=dir('*_tread_*_new.mat');

tread_SI=[];
tread_sparsity=[];
tread_width=[];
for i=1:length(list)
    load(list(i).name);
    
    tread_SI=[tread_SI analysis.SI];
    tread_sparsity=[tread_sparsity analysis.sparsity];
    tread_width=[tread_width analysis.width(analysis.pc_list)];
end

tab_width=cell2mat(tab_width');
tab_width=tab_width(:,1);

tread_width=cell2mat(tread_width');
tread_width=tread_width(:,1);

%% PC proportions
list=dir('analysis_*_new.mat');
for i=1:length(list)
    load(list(i).name);
    
    disp([list(i).name '     ' num2str(length(analysis.pc_list)) '/' num2str(size(analysis.deconv,2)) '     ' num2str(length(analysis.pc_list)/size(analysis.deconv,2))]);
end

%%
figure; 
subplot(1,3,1); hold on
[f,x]=ecdf(tab_SI);
plot(x,f,'r');
[f,x]=ecdf(tread_SI);
plot(x,f,'k');
xlim([0 6]);
ylim([0 1]);
axis square

subplot(1,3,2); hold on
[f,x]=ecdf(tab_sparsity);
plot(x,f,'r');
[f,x]=ecdf(tread_sparsity);
plot(x,f,'k');
xlim([0 1]);
ylim([0 1]);
axis square

subplot(1,3,3); hold on
[f,x]=ecdf(tab_width./80.*100);
plot(x,f,'r');
[f,x]=ecdf(tread_width./80.*100);
plot(x,f,'k');
xlim([0 100]);
ylim([0 1]);
axis square

%%
figure; 
subplot(1,3,1); hold on
% boxplot([tab_SI tread_SI],[ones(1,length(tab_SI)) 2.*ones(1,length(tread_SI))],'boxstyle','filled','colors','rk','plotstyle','compact');
boxplot([tab_SI tread_SI],[ones(1,length(tab_SI)) 2.*ones(1,length(tread_SI))],'colors','rk','symbol','o');
ylim([0 6]);
sigstar({[1,2]},1e-4);

subplot(1,3,2); hold on
boxplot([tab_sparsity tread_sparsity],[ones(1,length(tab_sparsity)) 2.*ones(1,length(tread_sparsity))],'colors','rk','symbol','o');
ylim([0 1]);
sigstar({[1,2]},1e-4);

subplot(1,3,3); hold on
boxplot([tab_width; tread_width],[ones(1,length(tab_width)) 2.*ones(1,length(tread_width))]','colors','rk','symbol','o');
ylim([0 100]);
sigstar({[1,2]},1);