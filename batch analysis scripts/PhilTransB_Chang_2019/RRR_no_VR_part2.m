%% load swr_classes

idx1 = [];
idx2 = [];

for f = 1:length(list)
    
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        full = fullfile(root(i).folder, root(i).name);
        
        clear h
        clear lfp
        load(fullfile(full, 'lfp1.mat'));
        idx1 = [idx1; lfp.swr_class];
        
        clear h
        clear lfp
        load(fullfile(full, 'lfp3.mat'));
        idx2 = [idx2; lfp.swr_class];
        
    end
end

%% err_map & sig_map for swr classes
figure('name','err r1')
subplot(2,3,1);
lol=r1_err_map;
lol=lol(idx1==-1,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR inhibited')
h(1)=subplot(2,3,4);
errorbar(mean(lol), sem(lol));

subplot(2,3,2);
lol=r1_err_map;
lol=lol(idx1==0,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR no modul')
h(2)=subplot(2,3,5);
errorbar(mean(lol), sem(lol));

subplot(2,3,3);
lol=r1_err_map;
lol=lol(idx1==1,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR excited')
h(3)=subplot(2,3,6);
errorbar(mean(lol), sem(lol));

linkaxes(h, 'y')

%
figure('name','err r2')
subplot(2,3,1);
lol=r2_err_map;
lol=lol(idx2==-1,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR inhibited')
h(1)=subplot(2,3,4);
errorbar(mean(lol), sem(lol));

subplot(2,3,2);
lol=r2_err_map;
lol=lol(idx2==0,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR no modul')
h(2)=subplot(2,3,5);
errorbar(mean(lol), sem(lol));

subplot(2,3,3);
lol=r2_err_map;
lol=lol(idx2==1,:);
[~,idx]=min(lol,[],2);
[~,idx]=sort(idx);
imagesc(lol(idx,:));
caxis([0 80]);
colorbar
title('SWR excited')
h(3)=subplot(2,3,6);
errorbar(mean(lol), sem(lol));

linkaxes(h, 'y')


%% decoding error comparison between SWR classes and rest1 rest2
lol1=r1_err_map;
lol1=lol1(idx1==0,:);
h(1)=errorshade(mean(lol1,'omitnan'), sem(lol1), 'color','k')
lol2=r1_err_map;
lol2=lol2(idx1==-1,:);
errorshade(mean(lol2,'omitnan'), sem(lol2), 'h',h(1), 'color','r')
s=ttest2(lol1,lol2,'tail','right','vartype','unequal');
idx = mean(lol1) + sem(lol1) + 10;
plot(find(~~s), idx(~~s), '*')

h(2)=errorshade(mean(lol1,'omitnan'), sem(lol1), 'color','k')
lol2=r1_err_map;
lol2=lol2(idx1==1,:);
errorshade(mean(lol2,'omitnan'), sem(lol2), 'h',h(2), 'color','y')
s=ttest2(lol1,lol2,'tail','right','vartype','unequal');
idx = mean(lol1) + sem(lol1) + 10;
plot(find(~~s), idx(~~s), '*')

lol1=r2_err_map;
lol1=lol1(idx2==0,:);
h(3)=errorshade(mean(lol1,'omitnan'), sem(lol1), 'color','k')
lol2=r2_err_map;
lol2=lol2(idx2==-1,:);
errorshade(mean(lol2,'omitnan'), sem(lol2), 'h',h(3), 'color','r')
s=ttest2(lol1,lol2,'tail','right','vartype','unequal');
idx = mean(lol1) + sem(lol1) + 10;
plot(find(~~s), idx(~~s), '*')

h(4)=errorshade(mean(lol1,'omitnan'), sem(lol1), 'color','k')
lol2=r2_err_map;
lol2=lol2(idx2==1,:);
errorshade(mean(lol2,'omitnan'), sem(lol2), 'h',h(4), 'color','y')
s=ttest2(lol1,lol2,'tail','right','vartype','unequal');
idx = mean(lol1) + sem(lol1) + 10;
plot(find(~~s), idx(~~s), '*')


linkaxes(h,'y')


%% min-max sorted err_map with no significant removed
[~,idx]=sort(mean(r1_err_map,2));
idx = intersect(idx, find(sum(r1_sig_map < .05, 2) > 0), 'stable');
figure
imagesc(r1_err_map(idx,:))
colormap jet
caxis([0 80]);
p = get(gca, 'position');
colorbar
set(gca, 'position', p);
figure
[~,idx]=sort(mean(r2_err_map,2));
idx = intersect(idx, find(sum(r2_sig_map < .05, 2) > 0), 'stable');
imagesc(r2_err_map(idx,:))
colormap jet
caxis([0 80]);
p = get(gca, 'position');
colorbar
set(gca, 'position', p);
