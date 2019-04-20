load lfp3.mat
load analysis.mat

lfp.detect_sce; ass=lfp.ensemble;
ass.cluster; ass.detect_sce;
ass.set_ops('order','cluster')

% EV=[]; REV=[];
for c=1:length(ass.clust)
    load 2\Plane1\deconv.mat
    thres = noRun(analysis.behavior.unit_vel);
    idx = analysis.behavior.unit_vel < thres;
    deconv(idx,:) = [];
    order=get_order(analysis);
    deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
    Rexp = corr(fast_smooth(deconv,30));

    load 1\Plane1\deconv.mat
    deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
    Rpre = corr(fast_smooth(deconv,30));

    load 3\Plane1\deconv.mat
    deconv = deconv(:,intersect(order,ass.clust{c},'stable'));
    Rpost = corr(fast_smooth(deconv,30));

    Rpre = triu(Rpre,1); Rpre = Rpre(Rpre~=0);
    Rexp = triu(Rexp,1); Rexp = Rexp(Rexp~=0);
    Rpost = triu(Rpost,1); Rpost = Rpost(Rpost~=0);

    PrePost = corr(Rpre, Rpost);
    ExpPre = corr(Rexp, Rpre);
    ExpPost = corr(Rexp, Rpost);

    EV = [EV ( (ExpPost - ExpPre*PrePost) / (sqrt((1 - ExpPre^2) * (1 - PrePost^2))) )^2];

    REV = [REV ( (ExpPre - ExpPost*PrePost) / (sqrt((1 - ExpPost^2) * (1 - PrePost^2))) )^2];
end

%%
figure
boxplot([EV' REV'])