layer=1;

index=[3 4 2 1];
count=1;
for i=layer:2:8
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_09',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_09',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[5 6 7 8 10 9];
count=1;
for i=layer:2:12
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_13',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_13',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[11 8 12];
count=1;
for i=layer:2:6
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_14',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_14',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[13 14 16 15];
count=1;
for i=layer:2:8
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_15',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_15',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[4];
count=1;
for i=layer:2:2
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_16',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_16',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[9];
count=1;
for i=layer:2:2
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_22',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis1=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc036','2017_11_22',num2str(i),2);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis2=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    analysis=merge_planes(analysis1,analysis2);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end








index=[3 4 2 1 3];
count=1;
for i=1:2:9
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc037','2017_11_09',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[5 6 7 9 10 8];
count=1;
for i=1:2:11
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc037','2017_11_13',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[11 4 16];
count=1;
for i=1:2:5
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc037','2017_11_14',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end

index=[14 15];
count=1;
for i=1:2:3
i
    [behavior,deconv,tcs,maskNeurons,mimg]=load_data('rsc037','2017_11_15',num2str(i),1);
    [behavior,deconv]=convert_behavior(behavior,tcs,deconv);
    analysis=pc_batch_analysis(behavior,deconv,maskNeurons,mimg);
    [centroids{index(count)},pc_list{index(count)}]=get_centroids(analysis);
    SI{index(count)}=analysis.SI;
    count=count+1;
end