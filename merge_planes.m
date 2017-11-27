function analysis=merge_planes(analysis1,analysis2,plotFlag)
% merge data from two planes, taking care to get rid of the same cells
if nargin<3
    plotFlag=0;
end

dist_thres=30; %pixels
sp_thres=.6;

mimg1=analysis1.mimg;
mimg2=analysis2.mimg;
% first, register mean images to correct translational differences
[optimizer,metric]=imregconfig('monomodal');
trans=imregtform(mimg2,mimg1,'translation',optimizer,metric);
reg=imwarp(mimg2,trans,'OutputView',imref2d(size(mimg1)));

if plotFlag
    figure
    imshowpair(mimg1,reg,'scaling','joint');
    title('overlaid planes');
end

trans.T=round(trans.T);
analysis2.maskNeurons=imwarp(analysis2.maskNeurons,trans,'OutputView',imref2d(size(mimg1)));
% get centroids and compute distributions of distance between nearest
% neighbours and other neighbours
cent1=get_centroids(analysis1);
cent2=get_centroids(analysis2);

x1=zeros(1,range(range(cent1)));
y1=zeros(1,range(range(cent1)));
x2=zeros(1,range(range(cent2)));
y2=zeros(1,range(range(cent2)));
for i=1:max(max(cent1))
    [x1(i),y1(i)]=find(cent1==i);
end

for i=1:max(max(cent2))
    [x2(i),y2(i)]=find(cent2==i);
end

dist=sqrt((x1'-x2).^2+(y1'-y2).^2);
% dist=dist(dist<thres);
% fit data to lognormal + logistic functions
% logi=@(x,L,k,x0) (L./(1+exp(-k.*(x-x0))));
% pdf=@(x,mu,sig,L,k,x0) lognpdf(x,mu,sig) + logi(x,L,k,x0);
% 
% [N,edges]=histcounts(dist,50);
% [~,idx]=sort(N);
% idx=edges(idx(end-1:end));
% mu=5;
% L=max(N);
% sig=1;
% k=1;
% x0=thres/2;
% start=[mu sig L k x0];
% 
% estimates=mle(dist,'pdf',pdf,'start',start);
dist=dist<dist_thres;

for i=1:size(dist,1)
    for j=1:size(dist,2)
        if dist(i,j) % distance too close?
            if corr2(analysis1.maskNeurons==i,analysis2.maskNeurons==j)>sp_thres %spatial corr too high?
                if any(analysis1.pc_list==i) || any(analysis2.pc_list==j) % are they even place cells?
                    if corr(analysis1.raw_stack(:,i),analysis2.raw_stack(:,j))>0.6
                        dist(i,j)=false;
                    end
                end
            end
        end
    end
end

dist=any(dist);
deconv=[analysis1.deconv analysis2.deconv(:,~dist)];
dist=find(dist);
same_cells=zeros(size(analysis1.maskNeurons));
count=0;
for i=dist
    same_cells(analysis2.maskNeurons==i-count)=10;
    analysis2.maskNeurons(analysis2.maskNeurons==i-count)=0;
    analysis2.maskNeurons(analysis2.maskNeurons>i-count)=analysis2.maskNeurons(analysis2.maskNeurons>i-count)-1;
    count=count+1;
end
analysis2.maskNeurons(analysis2.maskNeurons>0)=analysis2.maskNeurons(analysis2.maskNeurons>0)+range(range(analysis1.maskNeurons));
overlap=logical(logical(analysis1.maskNeurons).*logical(analysis2.maskNeurons));
new_mask=analysis1.maskNeurons+analysis2.maskNeurons;
new_mask(overlap)=analysis1.maskNeurons(overlap);
new_mimg=analysis1.mimg+analysis2.mimg;
analysis=pc_batch_analysis(analysis1.behavior,deconv,new_mask,new_mimg);

if plotFlag
    figure
    mask=~~analysis.maskNeurons+same_cells;
    mask(mask>=10)=2;
    imagesc(mask);
    axis square
    colormap hot
    title('overlapping cells');
end