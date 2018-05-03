function analysis=cat_analysis(behavior,deconv,bins)

if nargin<3
    bins=50;
end

% signal=ca_filt(deconv); 
signal=deconv;

idx=find(behavior.object_frame);

stack=zeros(length(idx)-1,bins,size(signal,2));
for i=1:length(idx)-1
    chunk=signal(idx(i):idx(i+1)-1,:);
    edges=[1:size(chunk,1)/bins:size(chunk,1) size(chunk,1)];
    edges=round(edges);
    
    temp=zeros(bins,size(signal,2));
    for j=1:length(edges)-1
        temp(j,:)=sum(chunk(edges(j):edges(j+1),:));
    end
    stack(i,:,:)=reshape(temp,1,bins,[]);
end

analysis.stack=stack;
analysis.bins=bins;