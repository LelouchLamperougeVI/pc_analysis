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
%   'par', true (default)
%       use parallel processing to speed up (can benefit very long
%       recordings)
%
%   'pic', false (default)
%       test for alternative hypothesis of "path integrator cells"

behavior=varargin{1};
deconv=varargin{2};

[maskFlag,testFlag,parFlag,shuffles,bins,sd,sig]=parse_input(varargin);

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

vr_length=round(range(unit_pos));
sr=1/mean(diff(frame_ts));

thres=noRun(unit_vel);

[psth,raw_psth,raw_count,edges,raw_stack,stack,Pi,vel_stack]=getStack(bins,sd,vr_length,deconv,thres,unit_pos,unit_vel,frame_ts,trials);

%SI test
signal=cell2mat(raw_count');
sizes=cellfun(@(x) size(x,1),raw_count);
SI=get_si(raw_count,edges,Pi);
if testFlag==1 || testFlag==3
    SI=[SI;zeros(shuffles,length(SI))];
    if parFlag
        parfor i=1:shuffles
%             temp=mat_circshift(signal,randi(size(signal,1),1,size(signal,2)));
            temp=signal(randperm(size(signal,1)),:);
            temp=mat2cell(temp,sizes,size(temp,2));
%             perm=randperm(numel(raw_count(:,:,1)));
%             perm=reshape(perm,size(raw_count,1),size(raw_count,2));
%             perm=repmat(perm,1,1,size(raw_count,3));
%             perm=perm+reshape(0:numel(raw_count(:,:,1)):numel(raw_count)-1,1,1,[]);
%             temp=raw_count(perm);
%             temp1=raw_psth(perm);
%             perm=randperm(length(raw_count));
            SI(i+1,:)=get_si(temp,edges,Pi);
        end
    else
        for i=1:shuffles
%             perm=randperm(numel(raw_count(:,:,1)));
%             perm=reshape(perm,size(raw_count,1),size(raw_count,2));
%             perm=repmat(perm,1,1,size(raw_count,3));
%             perm=perm+reshape(0:numel(raw_count(:,:,1)):numel(raw_count)-1,1,1,[]);
%             temp=raw_count(perm);
%             temp1=raw_psth(perm);
            perm=randperm(length(raw_count));
            SI(i+1,:)=get_si(raw_count(perm),edges,Pi);
        end
    end

    pval=1-sum(SI(1,:)>SI(2:end,:))./shuffles;
%     pval=sum(SI(1,:)>SI(2:end,:))./shuffles;
    pc_list=find(pval<sig);
    SI=SI(1,pc_list);
end


%PC width
width=cell(1,size(raw_stack,2));
for i=1:size(raw_stack,2)
    [pc_width,pc_loc]=ricker_test(stack(:,i));
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
% baseline_thres=range(raw_stack).*.2+min(raw_stack);
% width_series=raw_stack>baseline_thres;
% width_series=width_series(:,pc_list);
% 
% for i=1:size(width_series,2)
%     temp=width_series(:,i)';
%     start=strfind(temp,[0 1]);
%     ending=strfind(temp,[1 0]);
%     if temp(1)==1
%         start=[1 start];
%     end
%     if temp(end)==1
%         ending=[ending length(temp)];
%     end
%     temp=ending-start;
%     width(i)=max(temp);
% end
% width=width.*vr_length./bins;


%sparsity
Pi=Pi./sum(Pi);
sparsity=sum(Pi.*mean(raw_psth,1),2).^2./sum(Pi.*mean(raw_psth,1).^2);
sparsity=reshape(sparsity,1,[]);
sparsity=sparsity(pc_list);


if maskFlag
    maskNeurons=varargin{maskFlag+1};
    mimg=varargin{maskFlag+2};
    analysis=v2struct(vr_length,sr,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,behavior,maskNeurons,mimg);
else
    analysis=v2struct(vr_length,sr,psth,raw_psth,raw_stack,Pi,vel_stack,stack,SI,pval,pc_list,sparsity,width,deconv,behavior);
end


function [maskFlag,testFlag,parFlag,shuffles,bins,sd,sig]=parse_input(inputs)
maskFlag=0;
testFlag=1;
parFlag=true;
shuffles=1000;
sd=4;
bins=50;
sig=0.05;

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
    end
    idx=idx+1;
end



