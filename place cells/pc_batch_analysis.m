function analysis=pc_batch_analysis(varargin)
% analysis = pc_batch_analysis(behavior, deconv, params);
%
% Parameters
%   'mask', maskNeurons, mimg
%       required for merging two planes
%
%   'test',
%       'si'     SI shuffle test
%       'ricker' new method that convolves tuning curve with a series of
%                ricker wavelets
%       'mixed' (default) two tests combined
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
original_deconv=deconv;

ops = parse_input(varargin);
shuffles= ops.shuffles;

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

vr_length=round(range(unit_pos));
fs=1/median(diff(frame_ts));

thres=noRun(unit_vel);

thres=(unit_vel>thres | unit_vel<-thres) & (trials(1) < frame_ts & trials(end) > frame_ts);
unit_vel=unit_vel(thres);
unit_pos=unit_pos(thres);
frame_ts=frame_ts(thres);

try
    deconv=ca_filt(deconv);
catch
    warning('Unable to filter deconv time-courses. This data likely came from the new Suite2P.');
end
deconv=deconv(thres,:);

[psth,raw_psth,raw_stack,mu_fr,Pi,stack,vel_stack]=getStack(ops.bins,ops.sd,vr_length,deconv,unit_pos,unit_vel,frame_ts,trials);

silent=sum(deconv, 1, 'omitnan')==0; %cell that don't firing during the running epochs

[SI, SI_marge]=get_si_skaggs(raw_stack,mu_fr,Pi);
if ops.testFlag==1 || ops.testFlag==3
    SI=[SI;zeros(shuffles,length(SI))];
    shuff_stack = zeros(ops.bins,size(SI,2), shuffles);
    h=waitbar(0,'permutation testing...');
    count=1;
    if ops.parFlag
        dq=parallel.pool.DataQueue;
        afterEach(dq,@updateBar);
        parfor i=1:shuffles
            %             temp=randperm(size(deconv,1),size(deconv,2));
            %             temp=mat_circshift(deconv,temp);
            temp=burst_shuffler(deconv);
            [~,~,shuff_raw_stack,shuff_mu,shuff_pi,shuff_stack(:,:,i)]=getStack(ops.bins,ops.sd,vr_length,temp,unit_pos,unit_vel,frame_ts,trials);
            SI(i+1,:)=get_si_skaggs(shuff_raw_stack,shuff_mu,shuff_pi);
            send(dq,i);
        end
        close(h);
    else
        for i=1:shuffles
            %             temp=randperm(size(deconv,1),size(deconv,2));
            %             temp=mat_circshift(deconv,temp);
            temp=burst_shuffler(deconv);
            [~,~,shuff_raw_stack,shuff_mu,shuff_pi,shuff_stack(:,:,i)]=getStack(ops.bins,ops.sd,vr_length,temp,unit_pos,unit_vel,frame_ts,trials);
            SI(i+1,:)=get_si_skaggs(shuff_raw_stack,shuff_mu,shuff_pi);
            updateBar;
        end
        close(h);
    end
    
    pval=1-sum(SI(1,:)>SI(2:end,:))./shuffles;
    pc_list=find(pval<ops.sig);
    %     SI=SI(1,pc_list);
    SI=SI(1,:);
end


%PC width
width=cell(1,size(raw_stack,2));
rick_rejects=cell(1,size(raw_stack,2));
for i=1:size(raw_stack,2)
    if silent(i)
        width{i}=[];
        continue;
    end
    %     [pc_width,pc_loc]=ricker_test(stack(:,i),raw_psth(:,:,i),frac_trials,mad,width_thres,io_ratio);
    %     [pc_width,pc_loc,~,rick_rejects{i}]=ricker_test(stack(:,i),psth{i},frac_trials,mad,width_thres,io_ratio,consecutive);
    [pc_width,pc_loc,~,rick_rejects{i}]=ricker_test(stack(:,i),raw_psth(:,:,i),ops.frac_trials,ops.mad,ops.width,ops.io_ratio,ops.consecutive);
    width{i}=[pc_width pc_loc];
end
if ops.testFlag==2 || ops.testFlag==3
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
% sparsity=sparsity(pc_list);


if ops.maskFlag
    maskNeurons=varargin{ops.maskFlag+1};
    mimg=varargin{ops.maskFlag+2};
    analysis=v2struct(vr_length,fs,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,original_deconv,behavior,silent,rick_rejects,maskNeurons,mimg, shuff_stack, SI_marge);
else
    analysis=v2struct(vr_length,fs,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,original_deconv,behavior,silent,rick_rejects, shuff_stack, SI_marge);
end

    function updateBar(~)
        waitbar(count/shuffles,h);
        count=count+1;
    end
end


function ops = parse_input(inputs)
ops.maskFlag=0;
ops.testFlag=3;
ops.parFlag=true;
ops.shuffles=1000;
ops.sd=4;
ops.bins=50;
ops.sig=0.05;
ops.frac_trials=1/3;
ops.mad=3;
ops.width=[.05 .8];
ops.io_ratio=2.5;
ops.consecutive=true;

idx=3;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'mask'
            ops.maskFlag=idx;
            idx=idx+1;
        case 'test'
            switch inputs{idx+1}
                case 'si'
                    ops.testFlag=1;
                case 'ricker'
                    ops.testFlag=2;
                case 'mixed'
                    ops.testFlag=3;
                otherwise
                    error('not a valid test');
            end
        case 'sig'
            ops.sig=inputs{idx+1};
        case 'shuffles'
            ops.shuffles=inputs{idx+1};
        case 'frac_trials'
            ops.frac_trials=inputs{idx+1};
        case 'mad'
            ops.mad=inputs{idx+1};
        case 'width'
            ops.width=inputs{idx+1};
        case 'io_ratio'
            ops.io_ratio=inputs{idx+1};
        case 'consecutive'
            ops.consecutive=inputs{idx+1};
        case 'bins'
            ops.bins=inputs{idx+1};
        case 'sd'
            ops.sd=inputs{idx+1};
        case 'par'
            ops.parFlag=inputs{idx+1};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx=idx+2;
end
end
