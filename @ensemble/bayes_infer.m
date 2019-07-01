function [decoded, P] = bayes_infer(obj,varargin)
% Use Bayesian decoding to infer position encoded during rest
% name, value   pairs:
%   bins        number of spatial bins to decode (default - whatever used for pc_batch_analysis)
%   sd          gaussian smoothing sigma in centimeters (default - pc_batch_analysis)
%   tau         size of time window for decoding is seconds (default - sampling interval)
%   trainer     deconv for training (default - analysis.original_dconv)
%   tester      resting state deconv (default - obj.twop.deconv)
% Outputs:
%   decoded     decoded position

ops = []; trainer = []; tester = []; behavior = []; %declare and initialize global vars (I don't care what you tell me dad! I'm gonna use global variables!)
parse_inputs;

unit_pos=behavior.unit_pos;
unit_vel=behavior.unit_vel;
frame_ts=behavior.frame_ts;
trials=behavior.trials;

tau = ops.tau * obj.twop.fs; % time window length

thres=noRun(unit_vel);
thres=(unit_vel>thres | unit_vel<-thres) & (trials(1) < frame_ts & trials(end) > frame_ts);

artifact = [true diff(unit_vel)~=0];
thres = logical(thres .* artifact);
artifact = [true diff(unit_pos)~=0];
thres = logical(thres .* artifact);

unit_vel=unit_vel(thres);
unit_pos=unit_pos(thres);
frame_ts=frame_ts(thres);

trainer=ca_filt(trainer);
trainer=trainer(thres,:);
trainer=fast_smooth(trainer,ops.sig);

% tester(any(isnan(tester), 2), :) = [];
% tester=ca_filt(tester);
tester=fast_smooth(tester,ops.sig);
tester = movmean(tester,round(tau),1);

trainer= (trainer-min(trainer)) ./ range(trainer);
tester = (tester-min(tester)) ./ range(tester);

[~,~,stack]=getStack(ops.bins,ops.sd,obj.analysis.vr_length,trainer,unit_pos,unit_vel,frame_ts,trials);

% decoded = zeros(1 + length(obj.clust), size(tester,1));
% P = zeros(ops.bins, size(tester,1), 1 + length(obj.clust));
% count = 1;
decode(1:length(obj.analysis.psth));
% for i = 1:length(obj.clust)
%     count = count + 1;
%     decode(obj.clust{i});
% end

% figure;
% h(1)=subplot(2,1,1); imagesc(tester(:,obj.order)');
% colormap jet
% % h(1)=subplot(2,1,1); imagesc(zscore(tester(:,obj.analysis.order))');
% h(2)=subplot(2,1,2); imagesc(P(:,:,1));
% linkaxes(h, 'x');


    function decode(cluster)
        pr = prod(stack(:,cluster).^permute(tester(:,cluster),[3 2 1]), 2) .* exp(-tau .* sum(stack(:,cluster),2));
%         pr = prod(stack(:,cluster).^permute(tester(:,randperm(length(cluster))),[3 2 1]), 2) .* exp(-tau .* sum(stack(:,cluster),2));
        pr = squeeze(pr);
        P = pr ./ sum(pr,1);
        [~,decoded] = max(pr);
%         P(:,:,count) = pr ./ sum(pr,1);
%         [~,decoded(count, :)] = max(pr);
    end

    function parse_inputs
        ops.bins = size(obj.analysis.stack, 1);
        ops.sd = 4;
        ops.tau = 1 / obj.twop.fs;
        ops.sig = 5;
        
        trainer = obj.analysis.original_deconv;
        tester = obj.twop.deconv;
        behavior = obj.analysis.behavior;
        
        count = 1;
        while count < length(varargin)
            switch lower(varargin{count})
                case 'bins'
                    ops.bins = varargin{count+1};
                case 'sd'
                    ops.sd = varargin{count+1};
                case 'tau'
                    ops.tau = varargin{count+1};
                case {'train', 'trainer', 'training'}
                    trainer = varargin{count+1};
                case {'test', 'tester', 'testing'}
                    tester = varargin{count+1};
                case {'behaviour', 'behavior'}
                    behavior = varargin{count+1};
                otherwise
                    error(['''' varargin{count} ''' is not a valid parameter']);
            end
            count = count+2;
        end
    end
end