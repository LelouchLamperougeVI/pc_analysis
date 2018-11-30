function extract_traces(obj)

cam=obj.cam;
mask=obj.nose_mask;
traces=zeros(1,obj.num_frames);
parfor ii=1:obj.num_frames
    frame=cam.read(ii);
    frame=rgb2gray(frame);
    traces(ii)=mean(frame(mask));
end
obj.traces=traces;
obj.original_traces=traces;