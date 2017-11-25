function plot_mask(analysis)

pcMask=zeros(size(analysis.maskNeurons));
for i=analysis.pc_list
    pcMask(analysis.maskNeurons==i)=2;
end

figure
imagesc(~~analysis.maskNeurons+pcMask);
axis square
colormap hot