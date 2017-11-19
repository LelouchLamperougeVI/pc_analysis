function run_sequences=get_run_sequence(behavior,tcs,deconv)
% construct sequences from running-related transients for sequence alignment

[behavior,deconv]=convert_behavior(behavior,tcs,deconv);

analysis=pc_batch_analysis(behavior,deconv);

v2struct(analysis);
v2struct(behavior);

deconv=deconv(:,analysis.pc_list);

run_sequences=zeros(length(trials_ts)-1,size(deconv,2));

for i=1:length(trials_ts)-1
    seq=max(deconv(trials_ts(i):trials_ts(i+1)-1,:));
    [~,run_sequences(i,:)]=sort(seq);
end