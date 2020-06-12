function behavior=plot_beh_2d(behavior)

x=behavior.x;
y=behavior.y;
ts=behavior.ts;
reward=behavior.reward;
punish=behavior.punish;
trial=behavior.trial;
cone=behavior.cone;
sphere=behavior.sphere;

try
    obj_probe=behavior.obj_probe;
catch
    behavior.probe=false;
end

if ~behavior.probe
    figure;
    hold on;
    plot(cone(:,1),cone(:,2),'g^','markerfacecolor','g','markersize',20);
    plot(sphere(:,1),sphere(:,2),'ro','markerfacecolor','r','markersize',20);
    
    xlim([-100 100]);
    ylim([0 140]);
    
    for i=randperm(length(trial)-1,30) %randomly plot 10 trials
%     for i=length(trial)-20:length(trial)-1
%     for i=1:length(trial)-1
        plot(x(trial(i)+1:trial(i+1)),y(trial(i)+1:trial(i+1)),'k');
        idx=reward(reward>trial(i) & reward<=trial(i+1));
        plot(x(idx),y(idx),'kd','markerfacecolor','k','markersize',10);
        idx=punish(punish>trial(i) & punish<=trial(i+1));
        plot(x(idx),y(idx),'kx','linewidth',4,'markerfacecolor','k','markersize',10);
    end
else
    k=9;
    for i=1:length(obj_probe)-1
        if ~mod(i-1,k)
            figure;
        end
        subplot(3,3,mod(i-1,k)+1);
        hold on;
        plot(cone{i}(:,1),cone{i}(:,2),'g^','markerfacecolor','g','markersize',20);
        plot(sphere{i}(:,1),sphere{i}(:,2),'ro','markerfacecolor','r','markersize',20);
        temp=trial(ts(trial)>obj_probe(i) & ts(trial)<obj_probe(i+1));
        for j=1:length(temp)-1
            plot(x(temp(j)+1:temp(j+1)),y(temp(j)+1:temp(j+1)),'k');
            idx=reward(reward>temp(j) & reward<=temp(j+1));
            plot(x(idx),y(idx),'kd','markerfacecolor','k','markersize',10);
            idx=punish(punish>temp(j) & punish<=temp(j+1));
            plot(x(idx),y(idx),'kx','linewidth',4,'markerfacecolor','k','markersize',10);
        end
    end
    
    occ=zeros(2); %[left_occupancy right_occupancy
                  % left_null      right_null]
    tot_occ=zeros(1,2); % [left right]
    rewPun=zeros(2); %[left_rew     right_rew
                     % left_pun     right_pun]
    numTrials=zeros(1,2);
                     
    l_rw_zone = [-30 70; -10 210];
    r_rw_zone = [10 70; 30 210];
    
    border=@(x,y,boundary) sum(x>=boundary(1,1) & x<=boundary(2,1) & y>=boundary(1,2) & y<=boundary(2,2));
    
    figure;
    if sum(cone{1}(:,1))<0
        alt=0;
        first=false;
    else
        alt=1;
        first=true;
    end
    first_l=true;
    first_r=true;
    for i=1:length(obj_probe)-1
        try
            if i>1
                if ~all(all(cone{i}==cone{i-1}))
                    alt=alt+1;
                end
            end
            subplot(1,2,mod(alt,2)+1);
            if ((alt==0 && ~first) || ( first && alt==2 )) && first_l
                hold on;
                plot(cone{i}(:,1),cone{i}(:,2),'g^','markerfacecolor','g','markersize',20);
                plot(sphere{i}(:,1),sphere{i}(:,2),'ro','markerfacecolor','r','markersize',20);
                first_l=false;
            end
            if alt==1 && first_r
                hold on;
                plot(cone{i}(:,1),cone{i}(:,2),'g^','markerfacecolor','g','markersize',20);
                plot(sphere{i}(:,1),sphere{i}(:,2),'ro','markerfacecolor','r','markersize',20);
                first_r=false;
            end
            temp=trial(ts(trial)>obj_probe(i) & ts(trial)<obj_probe(i+1));
            for j=1:length(temp)-1
                plot(x(temp(j)+1:temp(j+1)),y(temp(j)+1:temp(j+1)),'k');
                idx=reward(reward>temp(j) & reward<=temp(j+1));
                plot(x(idx),y(idx),'kd','markerfacecolor','k','markersize',10);
                idx=punish(punish>temp(j) & punish<=temp(j+1));
                plot(x(idx),y(idx),'kx','linewidth',4,'markerfacecolor','k','markersize',10);
                
                tot_occ(mod(alt,2)+1)=temp(j+1)-temp(j)+1 + tot_occ(mod(alt,2)+1);
                occ(mod(alt,2)+1,1)=border(x(temp(j)+1:temp(j+1)), y(temp(j)+1:temp(j+1)), l_rw_zone) + occ(mod(alt,2)+1,1);
                occ(mod(alt+1,2)+1,2)=border(x(temp(j)+1:temp(j+1)), y(temp(j)+1:temp(j+1)), r_rw_zone) + occ(mod(alt+1,2)+1,2);
                
                if ~mod(alt,2)
                    idx=reward(reward>temp(j) & reward<=temp(j+1));
                else
                    idx=punish(punish>temp(j) & punish<=temp(j+1));
                end
                rewPun(mod(alt,2)+1,1)=border(x(idx), y(idx), l_rw_zone) + rewPun(mod(alt,2)+1,1);
                if mod(alt,2)
                    idx=reward(reward>temp(j) & reward<=temp(j+1));
                else
                    idx=punish(punish>temp(j) & punish<=temp(j+1));
                end
                rewPun(mod(alt+1,2)+1,2)=border(x(idx), y(idx), r_rw_zone) + rewPun(mod(alt+1,2)+1,2);
                
                numTrials(mod(alt,2)+1)=1 + numTrials(mod(alt,2)+1);
            end
        catch
        end
    end

    % behavior.occ=diff(occ([2 1],:))./tot_occ;
    behavior.occ=occ./tot_occ;    
    behavior.rewPun=rewPun./numTrials;
end

figure;
rew=zeros(ceil(length(trial)/10),1);
pun=rew;
count=1;
for i=1:10:length(trial)-10
    rew(count)=sum(reward>trial(i) & reward<=trial(i+10))/10;
    pun(count)=sum(punish>trial(i) & punish<=trial(i+10))/10;
    count=count+1;
end
rew(count)=sum(reward>trial(i+10) & reward<=trial(end))/10;
pun(count)=sum(punish>trial(i+10) & punish<=trial(end))/10;

plot([rew pun rew./pun]);