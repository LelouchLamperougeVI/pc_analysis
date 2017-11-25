function analysis=merge_planes(analysis1,analysis2)

mimg1=analysis1.mimg;
mimg2=analysis2.mimg;

mask2=analysis2.maskNeurons;

[optimizer,metric]=imregconfig('monomodal');
trans=imregtform(mimg2,mimg1,'translation',optimizer,metric);
reg=imwarp(mimg2,trans,'OutputView',imref2d(size(mimg1)));

figure
imshowpair(mimg1,reg,'scaling','joint');

mask2=imwarp(mask2,trans,'OutputView',imref2d(size(mimg1)));

cent1=get_centroids(analysis1);
cent2=get_centroids(analysis2,mask2);