function [decoded, P, pos, err] = bayes_infer(obj,varargin)
% Use Bayesian decoding to infer position encoded during rest
% name, value   pairs:
%   bins        number of spatial bins to decode (default - whatever used for pc_batch_analysis)
%   sd          gaussian smoothing sigma in centimeters (default - pc_batch_analysis)
%   tau         size of time window for decoding is seconds (default - sampling interval)
%   trainer     deconv for training (default - analysis.original_dconv)
%   tester      resting state deconv (default - obj.twop.deconv)
%   shuffle     whether to use `place field shuffle' - (default false)
%   validate    validation flag to test decoding performence on RUN data with even-odd trials - (default false; overrides trainer/tester)
%   circ        whether to use circular error estimate - (default false)
%   plotFlag
% Outputs:
%   decoded     decoded position
%   P           probability matrix
%   pos         actual animal position
%   err         decoding error and SEM in centimetres

ops = []; trainer = []; tester = []; behavior = []; pos = []; err = []; %declare and initialize global vars (I don't care what you tell me dad! I'm gonna use global variables!)
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
trainer= (trainer-min(trainer)) ./ range(trainer); % normalize between 0 and 1
%DO NOT zscore normalize, since negative firing rates are not acceptable

% tester(any(isnan(tester), 2), :) = [];
% tester=ca_filt(tester);
tester=fast_smooth(tester,ops.sig);
tester = movmean(tester,round(tau),1);
tester = (tester-min(tester)) ./ range(tester);

if ops.validate
    evenodd;
end

[~,~,stack]=getStack(ops.bins,ops.sd,obj.analysis.vr_length,trainer,unit_pos,unit_vel,frame_ts,trials);
if ops.shuffle
    stack = bcircshift(stack, randi(size(stack, 1), size(stack, 2), 1));
end

% decoded = zeros(1 + length(obj.ensembles.clust), size(tester,1));
% P = zeros(ops.bins, size(tester,1), 1 + length(obj.ensembles.clust));
% count = 1;
decode(1:length(obj.analysis.psth));
% for i = 1:length(obj.ensembles.clust)
%     count = count + 1;
%     decode(obj.ensembles.clust{i});
% end

if ops.validate
    if ops.circ
        err = [arrayfun(@(x) mean(min(...
                    [abs(decoded(pos == x) - pos(pos == x)) ; ...
                    ops.bins - abs(decoded(pos == x) - pos(pos == x))] ...
                    )), 1:ops.bins)
               arrayfun(@(x) sem(min(...
                    [abs(decoded(pos == x) - pos(pos == x)) ; ...
                    ops.bins - abs(decoded(pos == x) - pos(pos == x))] ...
                    )'), 1:ops.bins)
              ];
    else
        err = [arrayfun(@(x) mean(abs(decoded(pos==x) - pos(pos==x))), 1:ops.bins)
               arrayfun(@(x) sem(abs(decoded(pos==x) - pos(pos==x))'), 1:ops.bins)];
    end
    err = err' .* obj.analysis.vr_length ./ ops.bins;
end

if ops.plotFlag
    figure;
    h(1)=subplot(2,1,1); imagesc('xdata', obj.twop.ts, 'cdata',tester(:,obj.intern.order)');
    rbmap(h(1), 'cmap',hot, 'caxis', [0 max(tester(:))]);
    ylim(h(1), [1 length(obj.intern.order)]);
    h(2)=subplot(2,1,2); imagesc('xdata', obj.twop.ts, 'cdata',P(:,:,1));
    rbmap(h(2), 'cmap',hot, 'caxis', [0 1]);
    ylim(h(2), [1 size(P,1)]);
    linkaxes(h, 'x');
end


    function decode(cluster)
        pr = prod(stack(:,cluster) .^ permute(tester(:,cluster),[3 2 1]), 2) .* exp(-tau .* sum(stack(:,cluster),2));
        %         pr = prod(stack(:,cluster).^permute(tester(:,randperm(length(cluster))),[3 2 1]), 2) .* exp(-tau .* sum(stack(:,cluster),2));
        pr = squeeze(pr);
        %         P = pr;
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
        ops.plotFlag = false;
        ops.shuffle = false;
        ops.validate = false;
        ops.circ = false;
        
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
                case 'shuffle'
                    ops.shuffle = varargin{count+1};
                case {'train', 'trainer', 'training'}
                    trainer = varargin{count+1};
                case {'test', 'tester', 'testing'}
                    tester = varargin{count+1};
                case {'behaviour', 'behavior'}
                    behavior = varargin{count+1};
                case {'validate','validation'}
                    ops.validate = varargin{count+1};
                case {'circ', 'circular'}
                    ops.circ = varargin{count+1};
                case {'plot', 'plotflag'}
                    ops.plotFlag= varargin{count+1};
                otherwise
                    error(['''' varargin{count} ''' is not a valid parameter']);
            end
            count = count+2;
        end
    end

    function evenodd
        trials_idx = knnsearch(frame_ts', trials');
        
        even_trials = trials(2:2:end);
        odd_trials = trials(1:2:end);
        
        even_idx=[];
        odd_idx=[];
        for i=2:2:length(trials)-2
            even_idx = [even_idx trials_idx(i):trials_idx(i+1)];
        end
        try
            even_idx = [even_idx trials_idx(i+2):trials_idx(i+3)];
            even_trials = [even_trials trials(i+3)];
        catch
        end
        for i=1:2:length(trials)-2
            odd_idx = [odd_idx trials_idx(i):trials_idx(i+1)];
        end
        try
            odd_idx = [odd_idx trials_idx(i+2):trials_idx(i+3)];
            odd_trials = [odd_trials trials(i+3)];
        catch
        end
        
        tester = trainer(even_idx,:);
        pos = movmedian(unit_pos,round(tau),1);
        pos = pos(even_idx);
        pos = discretize(pos, linspace(min(pos), max(pos), ops.bins+1));
        
        trainer = trainer(odd_idx,:);
        unit_pos = unit_pos(odd_idx);
        unit_vel = unit_vel(odd_idx);
        frame_ts = frame_ts(odd_idx);
        trials = odd_trials;
    end
end