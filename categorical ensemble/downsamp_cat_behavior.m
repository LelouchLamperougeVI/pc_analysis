function behavior=downsamp_cat_behavior(behavior,bins)

idx=1:bins:length(behavior.frame_ts);
temp=knnsearch(idx',find(behavior.object_frame)');
behavior.object_frame=zeros(1,length(idx));
behavior.object_frame(temp)=1;
temp=knnsearch(idx',find(behavior.black_frame)');
behavior.black_frame=zeros(1,length(idx));
behavior.black_frame(temp)=1;
behavior.frame_ts=behavior.frame_ts(1:bins:end);
behavior.object_ts=behavior.frame_ts(logical(behavior.object_frame));
behavior.black_ts=behavior.frame_ts(logical(behavior.black_frame));