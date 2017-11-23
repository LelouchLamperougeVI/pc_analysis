function run_sequences=get_run_sequence(behavior,tcs,deconv)
% construct sequences from running-related transients for sequence alignment

% [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
% 
% analysis=pc_batch_analysis(behavior,deconv);
% 
% v2struct(analysis);
% v2struct(behavior);

% deconv=deconv(:,analysis.pc_list);

run_sequences=zeros(length(behavior.trials_ts)-1,size(deconv,2));

for i=1:length(behavior.trials_ts)-1
    seq=deconv(behavior.trials_ts(i):behavior.trials_ts(i+1)-1,:);
    [~,seq]=max(seq);
    [~,run_sequences(i,:)]=sort(seq);
end