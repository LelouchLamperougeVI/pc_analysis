function analysis=pc_batch_analysis(varargin)
% analysis = pc_batch_analysis(behavior, deconv, params);
%
% Parameters
%   'mask', maskNeurons, mimg
%       required for merging two planes
% 
%   'test',
%       'si' (default) SI shuffle test
%       'ricker' new method that convolves tuning curve with a series of
%               ricker wavelets
%       'mixed' two tests combined
%
%   'shuffles', 1000 (default)
%       number of shuffles for tests
%
%   'sig', 0.05 (default)
%       significance threshold for p-value
%
%   'bins', 50 (default)
%       number of spatial bins
%
%   'sd', 4 (default)
%       smoothing kernel s.d. in cm
%
%   'mad', 3 (default)
%       MAD threshold for Ricker test
%
%   'frac_trials', 1/3 (default)
%       fraction of active trials in place fields (for Ricker test)
%
%   'consecutive', true (default)
%       complementary parameter to frac_trials; whether trials need to be consecutive
%
%   'width', [.05 .8] (default)
%       minimum and maximum width of place fields expressed in fraction
%
%   'io_ratio', 2.5 (default)
%       minimum ratio between in-field vs out-field fr
%
%   'par', true (default)
%       use parallel processing to speed up (can benefit very long
%       recordings)
%
%   'pic', false (default)
%       test for alternative hypothesis of "path integrator cells"

behavior=varargin{1};
deconv=varargin{2};

[maskFlag,testFlag,parFlag,shuffles,bins,sd,sig,frac_trials,mad,width_thres,io_ratio,consecutive]=parse_input(varargin);

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

vr_length=round(range(unit_pos));
sr=1/mean(diff(frame_ts));

thres=noRun(unit_vel);

thres=(unit_vel>thres | unit_vel<-thres) & (trials(1) < frame_ts & trials(end) > frame_ts);
unit_vel=unit_vel(thres);
unit_pos=unit_pos(thres);
frame_ts=frame_ts(thres);

deconv=ca_filt(deconv);
deconv=deconv(thres,:);

[psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack]=getStack(bins,sd,vr_length,deconv,unit_pos,unit_vel,frame_ts,trials);

silent=sum(deconv)==0; %cell that don't firing during the running epochs

SI=get_si_skaggs(raw_stack,mu_fr,Pi);
if testFlag==1 || testFlag==3
    SI=[SI;zeros(shuffles,length(SI))];
    h=waitbar(0,'permutation testing that shit...');
    count=1;
    if parFlag
        dq=parallel.pool.DataQueue;
        afterEach(dq,@updateBar);
        parfor i=1:shuffles
%             temp=randperm(size(deconv,1),size(deconv,2));
%             temp=mat_circshift(deconv,temp);
            temp=burst_shuffler(deconv);
            [~,~,shuff_stack,shuff_mu,shuff_pi]=getStack(bins,sd,vr_length,temp,unit_pos,unit_vel,frame_ts,trials);
            SI(i+1,:)=get_si_skaggs(shuff_stack,shuff_mu,shuff_pi);
            send(dq,i);
        end
        close(h);
    else
        for i=1:shuffles
%             temp=randperm(size(deconv,1),size(deconv,2));
%             temp=mat_circshift(deconv,temp);
            temp=burst_shuffler(deconv);
            [~,~,shuff_stack,shuff_mu,shuff_pi]=getStack(bins,sd,vr_length,temp,unit_pos,unit_vel,frame_ts,trials);
            SI(i+1,:)=get_si_skaggs(shuff_stack,shuff_mu,shuff_pi);
            updateBar;
        end
        close(h);
    end

    pval=1-sum(SI(1,:)>SI(2:end,:))./shuffles;
    pc_list=find(pval<sig);
%     SI=SI(1,pc_list);
    SI=SI(1,:);
end


%PC width
width=cell(1,size(raw_stack,2));
for i=1:size(raw_stack,2)
    if silent(i)
        width{i}=[];
        continue;
    end
%     [pc_width,pc_loc]=ricker_test(stack(:,i),raw_psth(:,:,i),frac_trials,mad,width_thres,io_ratio);
    [pc_width,pc_loc]=ricker_test(stack(:,i),psth{i},frac_trials,mad,width_thres,io_ratio,consecutive);
    width{i}=[pc_width pc_loc];
end
if testFlag==2 || testFlag==3
    pc_list2=arrayfun(@(x) isempty(width{x}),1:size(raw_stack,2));
    pc_list2=find(~pc_list2);
    if exist('pc_list','var')
        pc_list=intersect(pc_list,pc_list2);
    else
        pc_list=pc_list2;
        pval=[];
    end
end

%sparsity
Pi=Pi./sum(Pi);
% sparsity=sum(Pi.*(mean(raw_psth,1).^2),2)./(mean(mean(raw_psth,1),2).^2);
sparsity=sum(Pi.*raw_stack).^2./sum(Pi.*raw_stack.^2);
% sparsity=shiftdim(sparsity,1);
sparsity=sparsity(pc_list);


if maskFlag
    maskNeurons=varargin{maskFlag+1};
    mimg=varargin{maskFlag+2};
    analysis=v2struct(vr_length,sr,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,behavior,silent,maskNeurons,mimg);
else
    analysis=v2struct(vr_length,sr,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,behavior,silent);
end

    function updateBar(~)
        waitbar(count/shuffles,h);
        count=count+1;
    end
end


function [maskFlag,testFlag,parFlag,shuffles,bins,sd,sig,frac_trials,mad,width,io_ratio,consecutive]=parse_input(inputs)
maskFlag=0;
testFlag=1;
parFlag=true;
shuffles=1000;
sd=4;
bins=50;
sig=0.05;
frac_trials=1/3;
mad=3;
width=[.05 .8];
io_ratio=2.5;
consecutive=true;

idx=3;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'mask'
            maskFlag=idx;
            idx=idx+2;
        case 'test'
            idx=idx+1;
            switch inputs{idx}
                case 'si'
                    testFlag=1;
                case 'ricker'
                    testFlag=2;
                case 'mixed'
                    testFlag=3;
                otherwise
                    error('not a valid test');
            end
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 'shuffles'
            idx=idx+1;
            shuffles=inputs{idx};
        case 'frac_trials'
            idx=idx+1;
            frac_trials=inputs{idx};
        case 'mad'
            idx=idx+1;
            mad=inputs{idx};
        case 'width'
            idx=idx+1;
            width=inputs{idx};
        case 'io_ratio'
            idx=idx+1;
            io_ratio=inputs{idx};
        case 'consecutive'
            idx=idx+1;
            consecutive=inputs{idx};
        case 'bins'
            idx=idx+1;
            bins=inputs{idx};
        case 'sd'
            idx=idx+1;
            sd=inputs{idx};
        case 'par'
            idx=idx+1;
            parFlag=inputs{idx};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx=idx+1;
end
end