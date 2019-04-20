%% RUN: remove non-running epochs
load 2\Plane1\deconv.mat
thres = noRun(analysis.behavior.unit_vel);
idx = analysis.behavior.unit_vel < thres;
deconv(idx,:) = [];
order=get_order(analysis);
% deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
deconv = deconv(:,intersect(order,cell2mat(ass.clust),'stable'));
Rexp = corr(fast_smooth(deconv,30));
figure
imagesc(Rexp); colormap jet

%%
load 1\Plane1\deconv.mat
% deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
deconv = deconv(:,intersect(order,cell2mat(ass.clust),'stable'));
Rpre = corr(fast_smooth(deconv,30));
figure
imagesc(Rpre); colormap jet

load 3\Plane1\deconv.mat
% deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
deconv = deconv(:,intersect(order,cell2mat(ass.clust),'stable'));
Rpost = corr(fast_smooth(deconv,30));
figure
imagesc(Rpost); colormap jet

%%
Rpre = triu(Rpre,1); Rpre = Rpre(Rpre~=0);
Rexp = triu(Rexp,1); Rexp = Rexp(Rexp~=0);
Rpost = triu(Rpost,1); Rpost = Rpost(Rpost~=0);

PrePost = corr(Rpre, Rpost);
ExpPre = corr(Rexp, Rpre);
ExpPost = corr(Rexp, Rpost);

EV = ( (ExpPost - ExpPre*PrePost) / (sqrt((1 - ExpPre^2) * (1 - PrePost^2))) )^2;

REV = ( (ExpPre - ExpPost*PrePost) / (sqrt((1 - ExpPost^2) * (1 - PrePost^2))) )^2;



