function extract_traces(obj)
% Extract traces of pixels intensity for individual masks

if ~isempty(obj.original_traces)
	warning('You already extracted the traces once. Are you sure you want to run it again? This may take a while. Do it in special circumstances, e.g. when you reselected masks.')
	uinp = input('Do you still want to extract? Y/N [Y] ', 's');
	if ~strcmpi(uinp, 'y') && ~isempty(uinp)
		return
	end
end

cam=obj.cam;
masks = obj.masks;
masks = reshape(masks, size(masks,1)*size(masks,2), 1, size(masks,3));
masks = squeeze(masks);
traces=zeros(obj.num_frames, length(obj.mask_labels));
trace_std=zeros(obj.num_frames, length(obj.mask_labels));
parfor ii=1:obj.num_frames
    frame=cam.read(ii);
    frame=double(rgb2gray(frame));
    traces(ii,:)=mean(frame(:) .* masks);
    trace_std(ii,:)=std(frame(:) .* masks);
end
obj.traces=traces;
obj.original_traces=traces;
obj.trace_std=trace_std;
