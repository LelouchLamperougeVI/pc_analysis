list = {'RSC036', 'RSC037', 'RSC038'};
bad = {'2017_08_11', '2017_08_18', '2017_09_18'};

sce1 = []; sce2 = [];
count = 1;
for f = 1:length(list)
    
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        if ~any(strcmp( fullfile(root(i).folder, root(i).name) , fullfile('E:\HaoRan\RRR\RSC038', bad)))
            full = fullfile(root(i).folder, root(i).name)
            
            clear lfp
            load(fullfile(full, 'lfp1.mat'));
            try sce1{count} = [lfp.clust_SCE.idx]; catch sce1{count} = []; end
            
            clear lfp
            load(fullfile(full, 'lfp3.mat'));
            try sce2{count} = [lfp.clust_SCE.idx]; catch sce2{count} = []; end
            
            count = count+1
        end
    end
end


%%
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 sum( P1(:, sce1{i}(:,j), i) > P1_95(:, sce1{i}(:,j), i) , 2)]; temp1 = temp1./sum(temp1); end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 sum( P2(:, sce2{i}(:,j), i) > P2_95(:, sce2{i}(:,j), i) , 2)]; temp2 = temp2./sum(temp2); end; end

h = errorshade(nanmean(temp1,2), sem(temp1, 2));
errorshade(nanmean(temp2,2), sem(temp2, 2), 'h',h, 'colour','r');

%%
% null = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); null = [null {P2_95(:, sce2{i}(:,j), i)}]; end; end
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 {P1(:, sce1{i}(:,j), i)}]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 {P2(:, sce2{i}(:,j), i)}]; end; end

%%
k = 10;
null = []; for j = 1:size(sce2{k}, 2); null = [null {P2_95(:, sce2{k}(:,j), k)}]; end
temp = []; for j = 1:size(sce2{k}, 2); temp = [temp {P2(:, sce2{k}(:,j), k)}]; end
for i = 1:length(null)
    h = errorshade( nanmean(null{i},2), sem(null{i},2) );
    errorshade( nanmean(temp{i},2), sem(temp{i},2), 'h',h, 'colour','r' );
    ylim([0 .2])
end

%%
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 P1(:, sce1{i}(:,j), i)]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 P2(:, sce2{i}(:,j), i)]; end; end
[~,idx]=min(temp2);
[~,idx]=sort(idx);
figure
imagesc(temp2(:,idx)')

%%
figure;
temp = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp = [temp P1(:, sce1{i}(:,j), i)]; end; end
hold on
plot(sum(temp'<.05) ./ sum(temp(:)<.05))
plot(sum(temp'<.01) ./ sum(temp(:)<.01))
plot(sum(temp'<.001) ./ sum(temp(:)<.001))
ylim([0 .12])

figure;
temp = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp = [temp P2(:, sce2{i}(:,j), i)]; end; end
hold on
plot(sum(temp'<.05) ./ sum(temp(:)<.05))
plot(sum(temp'<.01) ./ sum(temp(:)<.01))
plot(sum(temp'<.001) ./ sum(temp(:)<.001))
ylim([0 .12])

%%
temp2(temp2==0) = .001;
[r,p] = corr( repmat( distance_idx', size(temp2,2), 1 ), temp2(:))
[r,p] = corr( repmat( (49:-1:0)', size(stack,2), 1 ), stack(:))


%%
lag = 5;
lol = zeros(2*lag+1, size(temp,2), 4);
for i=1:4; lol(:,:,i) = cell2mat( arrayfun( @(x) xcorr(temp(:,x), belt_num==i, lag, 'unbiased'), 1:size(temp,2), 'uniformoutput',false) ); end

for i = 1:size(lol,2)
    idx = min(lol(:,i,:),[],1);
    [~,idx] = min(idx);
    lol(:,i,1) = lol(:,i, idx );
end
[~,idx]=min(temp);
[~,idx]=sort(idx);
figure
imagesc(lol(:,idx)')


%%
temp1 = cell(1, length(sce1));
temp2 = cell(1, length(sce2));
for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1{i} = [temp1{i} P1(:, sce1{i}(:,j), i)]; end; end
for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2{i} = [temp2{i} P2(:, sce2{i}(:,j), i)]; end; end
temp1 = cellfun( @(x) nanmean(x,2), temp1, 'uniformoutput',false );
temp2 = cellfun( @(x) nanmean(x,2), temp2, 'uniformoutput',false );
temp1( cellfun(@isempty, temp1) ) = [];
temp2( cellfun(@isempty, temp2) ) = [];
temp1 = cell2mat(temp1);
temp2 = cell2mat(temp2);

%% accumarray alpha levels vs distance from cue
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 P1(:, sce1{i}(:,j), i)]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 P2(:, sce2{i}(:,j), i)]; end; end

edges = [ logspace( log10(.05), log10(.001), 20) 0 ]; %edge bins of alpha values in logarithmic scale
idx = repmat( distance_idx', size(temp1,2), 1 );
idx2 = discretize(temp1(:), edges(end:-1:1));
idx(isnan(idx2))=[];
idx2(isnan(idx2))=[];

m1 = accumarray( [idx2 idx+1], 1);
% m1 = m1 ./ sum(m1,2);
m1 = m1 ./ accumarray(distance_idx'+1,1)';
m1 = m1 ./ sum(m1(:));
figure;
imagesc('xdata', min(idx):max(idx), 'ydata',edges(end:-1:1), 'cdata',m1);
axis square

edges = [ logspace( log10(.05), log10(.001), 20) 0 ];
idx = repmat( distance_idx', size(temp2,2), 1 );
idx2 = discretize(temp2(:), edges(end:-1:1));
idx(isnan(idx2))=[];
idx2(isnan(idx2))=[];

m2 = accumarray( [idx2 idx+1], 1);
% m2 = m2 ./ sum(m2,2);
m2 = m2 ./ accumarray(distance_idx'+1,1)';
m2 = m2 ./ sum(m2(:));
figure;
imagesc('xdata', min(idx):max(idx), 'ydata',edges(end:-1:1), 'cdata',m2);
axis square

clc
lregress(0:8, sum(m1) ./ sum(m1(:)));
xlim([-2 10]); ylim([.06 .16]);
[r,p] = corr((min(distance_idx):max(distance_idx))', (sum(m1) ./ sum(m1(:)))' )
lregress(0:8, sum(m2) ./ sum(m2(:)));
xlim([-2 10]); ylim([.06 .16]);
[r,p] = corr((min(distance_idx):max(distance_idx))', (sum(m2) ./ sum(m2(:)))' )


%% accumarray alpha levels vs distance from reward
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 P1(:, sce1{i}(:,j), i)]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 P2(:, sce2{i}(:,j), i)]; end; end

edges = [ logspace( log10(.05), log10(.001), 20) 0 ];
idx = repmat( (49:-1:0)', size(temp1,2), 1 );
idx2 = discretize(temp1(:), edges(end:-1:1));
idx(isnan(idx2))=[];
idx2(isnan(idx2))=[];

m1 = accumarray( [idx2 idx+1], 1);
figure;
imagesc('xdata', min(idx):max(idx), 'ydata',edges(end:-1:1), 'cdata',m1);
colormap bone
colorbar
axis square

edges = [ logspace( log10(.05), log10(.001), 20) 0 ];
idx = repmat( (49:-1:0)', size(temp2,2), 1 );
idx2 = discretize(temp2(:), edges(end:-1:1));
idx(isnan(idx2))=[];
idx2(isnan(idx2))=[];

m2 = accumarray( [idx2 idx+1], 1);
figure;
imagesc('xdata', min(idx):max(idx), 'ydata',edges(end:-1:1), 'cdata',m2);
colormap bone
colorbar
axis square

clc
figure
lregress((49:-1:0), sum(m1) ./ sum(m1(:)));
xlim([-10 60]); ylim([0 .05]);
[r,p] = corr((min(idx):max(idx))', (sum(m1) ./ sum(m1(:)))' )
lregress((49:-1:0), sum(m2) ./ sum(m2(:)));
xlim([-10 60]); ylim([0 .05]);
[r,p] = corr((min(idx):max(idx))', (sum(m2) ./ sum(m2(:)))' )


%%
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 {P1(:, sce1{i}(:,j), i)}]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 {P2(:, sce2{i}(:,j), i)}]; end; end

m1 = zeros(length(unique(distance_idx)), length(temp1));
for i = 1:length(temp1)
    edges = [ logspace( log10(.05), log10(.001), 10) 0 ];
    idx = repmat( distance_idx', size(temp1{i},2), 1 );
    idx2 = discretize(temp1{i}(:), edges(end:-1:1));
    idx(isnan(idx2))=[];
    idx2(isnan(idx2))=[];

    m = accumarray( [idx2 idx+1], 1, [length(edges)-1 length(unique(distance_idx))]);
    m = m ./ accumarray(distance_idx'+1,1)';
    m1(:,i) = sum(m) ./ sum(m(:));
end

m2 = zeros(length(unique(distance_idx)), length(temp2));
for i = 1:length(temp2)
    edges = [ logspace( log10(.05), log10(.001), 10) 0 ];
    idx = repmat( distance_idx', size(temp2{i},2), 1 );
    idx2 = discretize(temp2{i}(:), edges(end:-1:1));
    idx(isnan(idx2))=[];
    idx2(isnan(idx2))=[];

    m = accumarray( [idx2 idx+1], 1, [length(edges)-1 length(unique(distance_idx))]);
    m = m ./ accumarray(distance_idx'+1,1)';
    m2(:,i) = sum(m) ./ sum(m(:));
end


%%
temp1 = []; for i = 1:length(sce1); for j = 1:size(sce1{i}, 2); temp1 = [temp1 {P1(:, sce1{i}(:,j), i)}]; end; end
temp2 = []; for i = 1:length(sce2); for j = 1:size(sce2{i}, 2); temp2 = [temp2 {P2(:, sce2{i}(:,j), i)}]; end; end
m1 = cell2mat( cellfun(@(x) sum(x<.01,2) ./ size(x,2), temp1, 'uniformoutput',false) );
m2 = cell2mat( cellfun(@(x) sum(x<.01,2) ./ size(x,2), temp2, 'uniformoutput',false) );

errorshade( nanmean(m1,2), sem(m1,2))
errorshade( nanmean(m2,2), sem(m2,2))

[r,p] = corr( repmat( distance_idx', size(m2,2), 1 ), m2(:))


%% density place cells vs location (require width from part 2)
temp = cell2mat( cellfun(@(x) accumarray(x(:,2), 1/size(x,1), [50 1]), trials_width_clust1, 'uniformoutput',false) );
null = cell2mat( cellfun(@(x) accumarray(x(:,2), 1/size(x,1), [50 1]), trials_width_no_clust1, 'uniformoutput',false) );
idx = arrayfun(@(x) length(session(x).rest1.SI_clust_clust), 1:length(session));
null = cell2mat( arrayfun(@(x) repmat(null(:,x), 1, idx(x)), 1:length(idx), 'uniformoutput',false) );

errorshade( nanmean(temp-null,2), sem(temp-null,2) );
% h = zeros(1,size(temp,1)); for i = 1:size(temp,1); [~,h(i)] = ttest2(temp(i,:), 1/50, 'tail','right'); end


temp = cell2mat( cellfun(@(x) accumarray(x(:,2), 1/size(x,1), [50 1]), trials_width_clust2, 'uniformoutput',false) );
null = cell2mat( cellfun(@(x) accumarray(x(:,2), 1/size(x,1), [50 1]), trials_width_no_clust2, 'uniformoutput',false) );
idx = arrayfun(@(x) length(session(x).rest2.SI_clust_clust), 1:length(session));
null = cell2mat( arrayfun(@(x) repmat(null(:,x), 1, idx(x)), 1:length(idx), 'uniformoutput',false) );

errorshade( nanmean(temp-null,2), sem(temp-null,2) );


%%
temp1 = cell2mat(stack_clust2);
temp1 = [ max( cell2mat( arrayfun(@(x) xcorr(temp1(:,x), belt_num==1, 3, 'unbiased'), 1:size(temp1,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp1(:,x), belt_num==2, 3, 'unbiased'), 1:size(temp1,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp1(:,x), belt_num==3, 3, 'unbiased'), 1:size(temp1,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp1(:,x), belt_num==4, 3, 'unbiased'), 1:size(temp1,2), 'uniformoutput',false) ) ) ];
figure;
histogram(max(temp1),50,'normalization','probability');
% cdfplot(max(temp1));
hold on
temp2 = cell2mat(stack_clust1);
temp2 = [ max( cell2mat( arrayfun(@(x) xcorr(temp2(:,x), belt_num==1, 3, 'unbiased'), 1:size(temp2,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp2(:,x), belt_num==2, 3, 'unbiased'), 1:size(temp2,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp2(:,x), belt_num==3, 3, 'unbiased'), 1:size(temp2,2), 'uniformoutput',false) ) );
         max( cell2mat( arrayfun(@(x) xcorr(temp2(:,x), belt_num==4, 3, 'unbiased'), 1:size(temp2,2), 'uniformoutput',false) ) ) ];
histogram(max(temp2),50,'normalization','probability');
% cdfplot(max(temp2));




%% frac place cells
thres= 50/3; thres = 100;
temp1 = cell2mat( arrayfun(@(x) session(x).rest1.width, 1:length(session), 'uniformoutput',false)' );
temp1(temp1(:,1)>thres,:)=[];
temp1 = histcounts(temp1(:,2), 50);
temp2 = cell2mat( arrayfun(@(x) session(x).rest2.width, 1:length(session), 'uniformoutput',false)' );
temp2(temp2(:,1)>thres,:)=[];
temp2 = histcounts(temp2(:,2), 50);
temp3 = cell2mat( arrayfun(@(x) session(x).rest2.all_width, 1:length(session), 'uniformoutput',false)' );
temp3(temp3(:,1)>thres,:)=[];
temp3 = histcounts(temp3(:,2), 50);
figure
bar(temp1./sum(temp1)); title('R1');
% ylim([0 .25]);
figure
bar(temp2./sum(temp2)); title('R2');
% ylim([0 .25]);
figure
bar(temp3./sum(temp3)); title('pop');
% ylim([0 .25]);


temp1 = cell2mat( arrayfun(@(x) session(x).rest1.width_diff, 1:length(session), 'uniformoutput',false)' );
temp1(temp1(:,1)>thres,:)=[];
temp1 = histcounts(temp1(:,2), 50);
temp2 = cell2mat( arrayfun(@(x) session(x).rest2.width_diff, 1:length(session), 'uniformoutput',false)' );
temp2(temp2(:,1)>thres,:)=[];
temp2 = histcounts(temp2(:,2), 50);
temp3 = cell2mat( arrayfun(@(x) session(x).rest2.width_overlap, 1:length(session), 'uniformoutput',false)' );
temp3(temp3(:,1)>thres,:)=[];
temp3 = histcounts(temp3(:,2), 50);
figure
bar(temp1./sum(temp1)); title('lost');
% ylim([0 .25]);
figure
bar(temp2./sum(temp2)); title('gained');
% ylim([0 .25]);
figure
bar(temp3./sum(temp3)); title('stable');
% ylim([0 .25]);


%% SI
clc
figure
cdfplot(sum(SI_clust1))
hold on
cdfplot(sum(SI_clust2))
legend('r1','r2')
[h,p, stat] = kstest2( sum(SI_clust1), sum(SI_clust2) )
ranksum( sum(SI_clust1), sum(SI_clust2) )

figure
cdfplot(sum(SI_lost))
hold on
cdfplot(sum(SI_gained))
legend('lost','gained')
[h,p, stat] = kstest2( sum(SI_lost), sum(SI_gained) )
ranksum( sum(SI_lost), sum(SI_gained) )

%% de-marginalized SI
lims = [0 .07];
errorshade( nanmean(SI_clust1,2), sem(SI_clust1,2) );
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_clust1(i,:), SI_clust2(i,:), 'tail','right'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'k*');
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_clust1(i,:), SI_clust2(i,:), 'tail','left'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'r*');
ylim(lims);
title('r1');
errorshade( nanmean(SI_clust2,2), sem(SI_clust2,2) );
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_clust1(i,:), SI_clust2(i,:), 'tail','right'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'k*');
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_clust1(i,:), SI_clust2(i,:), 'tail','left'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'r*');
ylim(lims);
title('r2');
% errorshade( nanmean(SI_lost(:,sum(SI_lost)>1),2), sem(SI_lost(:,sum(SI_lost)>1),2) );
errorshade( nanmean(SI_lost,2), sem(SI_lost,2) );
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_lost(i,:), SI_gained(i,:), 'tail','right'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'k*');
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_lost(i,:), SI_gained(i,:), 'tail','left'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'r*');
ylim(lims);
title('lost');
% errorshade( nanmean(SI_gained(:,sum(SI_gained)>1),2), sem(SI_gained(:,sum(SI_gained)>1),2) );
errorshade( nanmean(SI_gained,2), sem(SI_gained,2) );
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_lost(i,:), SI_gained(i,:), 'tail','right'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'k*');
h=[]; for i = 1:size(SI_clust2,1); h(i) = ranksum(SI_lost(i,:), SI_gained(i,:), 'tail','left'); end
plot(find(h<.05), zeros(1,sum(h<.05)), 'r*');
ylim(lims);
title('gained');
errorshade( nanmean(SI_stable,2), sem(SI_stable,2) );
ylim(lims);
title('stable');



%% inter event interval histograms
sce_isi1(isinf(log(sce_isi1)))=[];
sce_isi2(isinf(log(sce_isi2)))=[];

edges = logspace(0,4,30);
d = histcounts(sce_isi1, edges);
figure
bar(edges(2:end), d./sum(d));
ylim([0 .2])
set(gca,'xscale','log')

d = histcounts(sce_isi2, edges);
figure
bar(edges(2:end), d./sum(d));
ylim([0 .2])
set(gca,'xscale','log')


edges = linspace(0,500,30);
d = histcounts(sce_isi1, edges);
figure
bar(edges(2:end), d./sum(d));
ylim([0 .2])

d = histcounts(sce_isi2, edges);
figure
bar(edges(2:end), d./sum(d));
ylim([0 .2])

