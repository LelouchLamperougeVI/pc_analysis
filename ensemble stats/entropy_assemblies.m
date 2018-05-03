function assemblies=entropy_assemblies(deconv,max_win,analysis)
% max_win: max number of frames for activation
if nargin<3
    ordered=1:size(deconv,2);
else
    bins=length(analysis.Pi);
    stack=analysis.raw_stack;

    stack=(stack-repmat(min(stack),bins,1));
    stack=stack./repmat(max(stack),bins,1);

    [~,idx]=max(stack);
    [~,ordered]=sort(idx);
end

deconv_mask=ca_filt(deconv)>0;
deconv=zscore(deconv);

scale=floor(log2(max_win));
sce=entropy_sce(deconv,scale,1);

series=zeros(length(find(sce)),size(deconv,2));
window_mask=false(size(deconv,1),length(find(sce)));

count=1;
for i=1:size(sce,2)
    for j=find(sce(:,i))'
        idx=j-floor(2^(i-1)/2):j+floor(2^(i-1)/2);
        series(count,:)=(mean(deconv(idx,:),1));
        window_mask(idx,count)=true;
        count=count+1;
    end
end

[assemblies,R]=lopes_pca(series,0,1);

% activation=false(size(deconv));
figure;
imagesc(fast_smooth(deconv(:,ordered),5)');
colormap gray
hold on
for i=1:length(assemblies)
    activation_mask=R(:,i)>prctile(R(:,i),95);
    activation_mask=window_mask(:,activation_mask);
    activation_mask=logical(sum(activation_mask,2));
    
    assembly_mask=false(1,size(deconv,2));
    assembly_mask(assemblies{i})=true;
    
    activation=activation_mask&assembly_mask&deconv_mask;
    [x,y]=find(activation(:,ordered));
    plot(x,y,'.');
end