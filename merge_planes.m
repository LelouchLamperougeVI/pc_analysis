function analysis=merge_planes(analysis1,analysis2)

mimg1=analysis1.mimg;
mimg2=analysis2.mimg;

[optimizer,metric]=imregconfig('monomodal');
trans=imregtform(mimg2,mimg1,'translation',optimizer,metric);
reg=imwarp(mimg2,trans,'OutputView',imref2d(size(mimg1)));

figure
imshowpair(mimg1,reg,'scaling','joint');

analysis2.maskNeurons=imwarp(analysis2.maskNeurons,trans,'OutputView',imref2d(size(mimg1)));

cent1=get_centroids(analysis1);
cent2=get_centroids(analysis2);

