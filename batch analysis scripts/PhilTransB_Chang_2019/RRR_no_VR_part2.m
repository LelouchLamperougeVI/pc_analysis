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


%% overlap + similarity
% ensembles have different degrees of overlap between r1 and r2
% which ensemble corresponds with which and to what degree are they
% overlapping?

r1_overlap = [];
r2_overlap = [];
for f = 1:length(list)
    
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        full = fullfile(root(i).folder, root(i).name)
        
        clear lfp
        load(fullfile(full, 'lfp1.mat'));
        ass1 = lfp.clust;
        
        clear lfp
        load(fullfile(full, 'lfp3.mat'));
        ass3 = lfp.clust;
        
        for ii = 1:length(ass1)
            idx = cellfun( @(x) length( intersect(ass1{ii}, x) ) ./ ( length(x) + length(ass1{ii}) - length( intersect(ass1{ii}, x) ) ), ass3 );
            [l, idx] = max(idx);
            r1_overlap = [r1_overlap; [l idx] ];
        end
        
        for ii = 1:length(ass3)
            idx = cellfun( @(x) length( intersect(ass3{ii}, x) ) ./ ( length(x) + length(ass3{ii}) - length( intersect(ass3{ii}, x) ) ), ass1 );
            [l, idx] = max(idx);
            r2_overlap = [r2_overlap; [l idx] ];
        end
    end
end


%% Pop stack mean
figure

h(1) = subplot(2,2,1);
stack=arrayfun(@(x) session(x).rest1.stack_clust_pc, 1:length(session), 'uniformoutput',false);
% stack=arrayfun(@(x) session(x).rest1.stack_clust_pc_z, 1:length(session), 'uniformoutput',false);
stack = cell2mat(stack);
% stack = fast_smooth(stack, 4/148*50);
% stack = (stack - min(stack)) ./ range(stack);
% stack = zscore(stack);
errorshade( mean(stack,2), sem(stack,2), 'h',h(1) );
title('R1 Clust PC');

h(2) = subplot(2,2,2);
stack=arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false);
% stack=arrayfun(@(x) session(x).rest2.stack_clust_pc_z, 1:length(session), 'uniformoutput',false);
stack = cell2mat(stack);
% stack = fast_smooth(stack, 4/148*50);
% stack = (stack - min(stack)) ./ range(stack);
% stack = zscore(stack);
errorshade( mean(stack,2), sem(stack,2), 'h',h(2) );
title('R2 Clust PC');

h(3) = subplot(2,2,3);
stack=arrayfun(@(x) session(x).rest1.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
% stack=arrayfun(@(x) session(x).rest1.stack_no_clust_pc_z, 1:length(session), 'uniformoutput',false);
stack = cell2mat(stack);
% stack = fast_smooth(stack, 4/148*50);
% stack = (stack - min(stack)) ./ range(stack);
% stack = zscore(stack);
errorshade( mean(stack,2), sem(stack,2), 'h',h(3) );
title('R1 No Clust PC');

h(4) = subplot(2,2,4);
stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
% stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc_z, 1:length(session), 'uniformoutput',false);
stack = cell2mat(stack);
% stack = fast_smooth(stack, 4/148*50);
% stack = (stack - min(stack)) ./ range(stack);
% stack = zscore(stack);
errorshade( mean(stack,2), sem(stack,2), 'h',h(4) );
title('R2 No Clust PC');

linkaxes(h, 'y')


%% pc width
stack_shuff = [];
width = [];
SI = [];
SI_all = [];
SI_clust1 = [];
SI_clust2 = [];
SI_no_clust1 = [];
SI_no_clust2 = [];
trials_stack = {};
trials_width_clust1 = {};
trials_width_no_clust1 = {};
trials_width_clust2 = {};
trials_width_no_clust2 = {};
vel_stack=[];
% bad = {};
bad = {'2017_08_11', '2017_08_18', '2017_09_18'};
% bad = {'2017_08_11', '2017_08_18', '2017_09_18', '2017_11_02'};
% bad = {'2017_08_11', '2017_08_18', '2017_09_13', '2017_09_18', '2017_11_02'};
count = 0;
for f = 1:length(list)
    
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        if ~any(strcmp( fullfile(root(i).folder, root(i).name) , fullfile('E:\HaoRan\RRR\RSC038', bad)))
            full = fullfile(root(i).folder, root(i).name)

            clear analysis
            load(fullfile(full, 'analysis.mat'));
            width = [width analysis.width(analysis.pc_list)];
            temp = analysis.SI(analysis.pc_list);
            idx = [cellfun(@(x) any( ismember( x(:,2), find(belt_num == 1) ) ), analysis.width(analysis.pc_list));...
                cellfun(@(x) any( ismember( x(:,2), find(belt_num == 2) ) ), analysis.width(analysis.pc_list));...
                cellfun(@(x) any( ismember( x(:,2), find(belt_num == 3) ) ), analysis.width(analysis.pc_list));...
                cellfun(@(x) any( ismember( x(:,2), find(belt_num == 4) ) ), analysis.width(analysis.pc_list)) ];
            SI = [SI, temp.*idx];
            idx = cellfun(@(x) arrayfun(@(y) any(ismember( x(:,2), y )), 1:50 ), analysis.width(analysis.pc_list), 'uniformoutput',false);
            idx = cell2mat(idx')';
            SI_all = [SI_all temp.*idx];
            trials_stack = [trials_stack {analysis.stack(:,analysis.pc_list)}];
            vel_stack = [vel_stack mean(analysis.vel_stack)'];
            stack_shuff = [stack_shuff mean( analysis.shuff_stack(:,analysis.pc_list,:), 3)];
            
            clear lfp
            load(fullfile(full, 'lfp1.mat'));
            trials_width_clust1 = [trials_width_clust1 cell2mat( analysis.width( intersect(analysis.pc_list, cell2mat(lfp.clust) ) )' )];
            trials_width_no_clust1 = [trials_width_no_clust1 cell2mat( analysis.width( intersect(analysis.pc_list, setxor( 1:length(analysis.psth), cell2mat(lfp.clust) ) ) )' )];
            temp = analysis.SI_marge(:, intersect(analysis.pc_list, cell2mat(lfp.clust) ) );
            SI_clust1 = [SI_clust1, temp];
            temp = analysis.SI_marge(:, intersect(analysis.pc_list, setxor( 1:length(analysis.psth), cell2mat(lfp.clust) ) ) );
            SI_no_clust1 = [SI_no_clust1, temp];
            
            clear lfp
            load(fullfile(full, 'lfp3.mat'));
            trials_width_clust2 = [trials_width_clust2 cell2mat( analysis.width( intersect(analysis.pc_list, cell2mat(lfp.clust) ) )' )];
            trials_width_no_clust2 = [trials_width_no_clust2 cell2mat( analysis.width( intersect(analysis.pc_list, setxor( 1:length(analysis.psth), cell2mat(lfp.clust) ) ) )' )];
            temp = analysis.SI_marge(:, intersect(analysis.pc_list, cell2mat(lfp.clust) ) );
            SI_clust2 = [SI_clust2, temp];
            temp = analysis.SI_marge(:, intersect(analysis.pc_list, setxor( 1:length(analysis.psth), cell2mat(lfp.clust) ) ) );
            SI_no_clust2 = [SI_no_clust2, temp];
            
            count = count+1
        end
    end
end

width = cell2mat(width');
trials_width_clust1 = cell2mat(trials_width_clust1');
trials_width_no_clust1 = cell2mat(trials_width_no_clust1');
trials_width_clust2 = cell2mat(trials_width_clust2');
trials_width_no_clust2 = cell2mat(trials_width_no_clust2');
SI(~SI) = nan;
SI_all(~SI_all) = nan;
SI_clust1(~SI_clust1) = nan;
SI_clust2(~SI_clust2) = nan;
SI_no_clust1(~SI_no_clust1) = nan;
SI_no_clust2(~SI_no_clust2) = nan;
%% cont'd
figure
h(1) = subplot(2,2,2);
p = accumarray(width, 1 / size(width,1) );
p = p(min(width(:,1)):2:end, :);
p = fast_smooth2(p, 'method','gauss', 'sig',1);
p = p ./ sum(p(:));
scale = unique(width(:,1));
scale = scale .* 3;
x = linspace(0, 150, range(width(:,2)) + 1);
imagesc('xdata',x, 'ydata', scale, 'cdata', p);
caxis([0 .02]);
axis equal
pos = get(h(1), 'position');
colorbar
cmap = bone; colormap(cmap(end:-1:1, :));
set( h(1), 'position',pos);

h(2) = subplot(2,2,1); %marginal y
plot(sum(p,2), scale);
pbaspect(h(2), pbaspect(h(1)) );

h(3) = subplot(2,2,4); %marginal x
plot(x, sum(p));
% axis square

linkaxes(h(1:2), 'y');
linkaxes(h([1 3]), 'x');

subplot(h(3));
xlim([0 150]); ylim([0 .05]);
subplot(h(2));
ylim([scale(1) scale(end)]); xlim([0 .2]);


%% cont'd all conditions
figure; hold on;

p = accumarray(trials_width_clust1, 1 / size(trials_width_clust1,1) );
p = p(min(trials_width_clust1(:,1)):2:end, :);
p = fast_smooth2(p, 'method','gauss', 'sig',1);
p = p ./ sum(p(:));
plot(sum(p));

p = accumarray(trials_width_no_clust1, 1 / size(trials_width_no_clust1,1) );
p = p(min(trials_width_no_clust1(:,1)):2:end, :);
p = fast_smooth2(p, 'method','gauss', 'sig',1);
p = p ./ sum(p(:));
plot(sum(p));

p = accumarray(trials_width_clust2, 1 / size(trials_width_clust2,1) );
p = p(min(trials_width_clust2(:,1)):2:end, :);
p = fast_smooth2(p, 'method','gauss', 'sig',1);
p = p ./ sum(p(:));
plot(sum(p));

p = accumarray(trials_width_no_clust2, 1 / size(trials_width_no_clust2,1) );
p = p(min(trials_width_no_clust2(:,1)):2:end, :);
p = fast_smooth2(p, 'method','gauss', 'sig',1);
p = p ./ sum(p(:));
plot(sum(p));

%% activity box plots
stack=[ arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false),...
        arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false) ];
stack=cell2mat(stack);
x = [];
g = [];
mu=[];err=[];
idx = [3 4 1 2];
for i = 1:4
    temp = stack(belt_num == idx(i), :);
    mu(i) = mean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];
end
i=5;
    temp = stack(1, :);
    mu(i) = mean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];
    
figure
errorbar(mu, err);
hold on
bar(mu);
ylim([.25 .31])

[~,~,stats] = anova1(x, g);
c = multcompare(stats);

%% activity error plots cont'd
x = [];
g = [];
mu=[];
err=[];
idx = [3 4 1 2];

stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false); stack=cell2mat(stack); stack = (stack - min(stack)) ./ range(stack);
for i = 1:4
    temp = stack(belt_num == idx(i), :);
    mu(i) = mean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];
    count = count+1;
end
i=5;
    temp = stack(1, :);
    mu(i) = mean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];

[~,~,stats] = anovan(x, g);
c = multcompare(stats)

figure
errorbar(mu, err);
hold on
bar(mu);
ylim([.23 .33])

%% clust x noclust correlogram
stack_clust = arrayfun(@(x) session(x).rest1.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack_no_clust = arrayfun(@(x) session(x).rest1.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);

stack_clust = cellfun(@(x) mean( x ,2), stack_clust, 'uniformoutput',false);
stack_no_clust = cellfun(@(x) mean( x ,2), stack_no_clust, 'uniformoutput',false);

stack_clust = cell2mat(stack_clust);
stack_no_clust = cell2mat(stack_no_clust);

[r,p] = corr(stack_clust', stack_no_clust');
p = p<.05;
p = bwboundaries(p);

figure('name','R1');
imagesc(r);
rbmap('caxis',[-.5 1], 'interp',120);
hold on
for i = 1:length(p)
    plot(p{i}(:,2), p{i}(:,1), 'k');
end


stack_clust = arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack_no_clust = arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);

stack_clust = cellfun(@(x) mean( x ,2), stack_clust, 'uniformoutput',false);
stack_no_clust = cellfun(@(x) mean( x ,2), stack_no_clust, 'uniformoutput',false);

stack_clust = cell2mat(stack_clust);
stack_no_clust = cell2mat(stack_no_clust);

[r,p] = corr(stack_clust', stack_no_clust');
p = p<.05;
p = bwboundaries(p);

figure('name','R2');
imagesc(r);
rbmap('caxis',[-.5 1], 'interp',120);
hold on
for i = 1:length(p)
    plot(p{i}(:,2), p{i}(:,1), 'k');
end



%% Cue xcorr
stack=arrayfun(@(x) session(x).rest1.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(stack);
r = cell2mat( arrayfun( @(x) xcorr(stack(:,x), (belt_num==4)', 5, 'coeff'), 1:size(stack,2), 'uniformoutput',false) );
errorshade(mean(r,2), sem(r,2))


%% Cue and/or reward correlation
clc
stack=arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false);
% stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
% stack = cell2mat(stack);
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:))
[r,p] = corr( repmat( distance_idx', size(stack,2), 1 ), stack(:))
% [r,p] = corr( repmat( distance_cues, size(stack,2), 1 ), stack(:))

%%
figure
clc
stack=arrayfun(@(x) session(x).rest1.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(1) = subplot(2,2,1);
lregress(repmat( (49:-1:0)', size(stack,2), 1 ), stack(:), 1, h(1));
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(2) = subplot(2,2,2);
lregress(repmat( (49:-1:0)', size(stack,2), 1 ), stack(:), 1, h(2));
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest1.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(3) = subplot(2,2,3);
lregress(repmat( (49:-1:0)', size(stack,2), 1 ), stack(:), 1, h(3));
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(4) = subplot(2,2,4);
lregress(repmat( (49:-1:0)', size(stack,2), 1 ), stack(:), 1, h(4));
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);

linkaxes(h, 'xy')

%%
figure
clc
stack=arrayfun(@(x) session(x).rest1.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(1) = subplot(2,2,1);
lregress(repmat( distance_idx', size(stack,2), 1 ), stack(:), 1, h(1));
[r,p] = corr( repmat( distance_idx', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest2.stack_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(2) = subplot(2,2,2);
lregress(repmat( distance_idx', size(stack,2), 1 ), stack(:), 1, h(2));
[r,p] = corr( repmat( distance_idx', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest1.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(3) = subplot(2,2,3);
lregress(repmat( distance_idx', size(stack,2), 1 ), stack(:), 1, h(3));
[r,p] = corr( repmat( distance_idx', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);
stack=arrayfun(@(x) session(x).rest2.stack_no_clust_pc, 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
h(4) = subplot(2,2,4);
lregress(repmat( distance_idx', size(stack,2), 1 ), stack(:), 1, h(4));
[r,p] = corr( repmat( distance_idx', size(stack,2), 1 ), stack(:));
disp(['r = ' num2str(r) ' p = ' num2str(p)]);

linkaxes(h, 'xy')


%% Speed x FR
% a sanity check...
stack = arrayfun(@(x) [session(x).rest2.stack_clust_pc session(x).rest2.stack_no_clust_pc], 1:length(session), 'uniformoutput',false);
stack = cell2mat(cellfun(@(x) mean(x,2), stack, 'uniformoutput',false));
[r,p] = corr(vel_stack(~isnan(vel_stack)), stack(~isnan(vel_stack)))
figure
plot(vel_stack(~isnan(vel_stack)), stack(~isnan(vel_stack)), '.k');
axis square


%% decoding error

errorshade(nanmean(r1_err_map)', sem(r1_err_map)');
errorshade(nanmean(r1_err_noclust)', sem(r1_err_noclust)')
errorshade(nanmean(r1_err_null)', sem(r1_err_null)')

errorshade(nanmean(r2_err_map)', sem(r2_err_map)');
errorshade(nanmean(r2_err_noclust)', sem(r2_err_noclust)')
errorshade(nanmean(r2_err_null)', sem(r2_err_null)')


%% SI
figure

h(1) = subplot(2,2,1);
errorshade( nanmean(SI_clust1,2), sem(SI_clust1,2), 'h',h(1) );
title('R1 Clust PC');

h(2) = subplot(2,2,2);
errorshade( nanmean(SI_clust2,2), sem(SI_clust2,2), 'h',h(2) );
title('R2 Clust PC');

h(3) = subplot(2,2,3);
errorshade( nanmean(SI_no_clust1,2), sem(SI_no_clust1,2), 'h',h(3) );
title('R1 No Clust PC');

h(4) = subplot(2,2,4);
errorshade( nanmean(SI_no_clust2,2), sem(SI_no_clust2,2), 'h',h(4) );
title('R2 No Clust PC');

linkaxes(h, 'y')


%% boxplots
tmp = SI_clust2;
x = [];
g = [];
mu=[];err=[];
idx = [3 4 1 2];
for i = 1:4
    temp = tmp(belt_num == idx(i), :);
    mu(i) = nanmean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];
end
i=5;
    temp = tmp(1, :);
    mu(i) = nanmean(temp(:));
    err(i) = sem(temp(:));
    x = [x; temp(:)];
    g = [g; i.* ones(numel(temp), 1)];
    
figure
errorbar(mu, err);
hold on
bar(mu);
ylim([0 .04])

[~,~,stats] = anova1(x, g);
c = multcompare(stats)















