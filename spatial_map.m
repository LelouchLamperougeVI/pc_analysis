bins=10;
winSize=800;
numWin=[6 3];
smoothing=2;

fov=835.76; %um
bregma=3; %top left corner

map=zeros(bins*numWin(1),bins*numWin(2));
SI_map=map;
mimg_map=zeros(winSize*numWin(1),winSize*numWin(2));

index=[13 12 0
    4 8 11
    3 7 10
    2 6 9
    1 5 14
    16 15 0];

mimg_index=repelem(index,winSize,winSize);
map_index=repelem(index,bins,bins);

for i=1:length(centroids)
    if ~isempty(centroids{i})
        all_cells=logical(centroids{i});
        place_cells=zeros(size(centroids{i}));
        SI_cells=zeros(size(centroids{i}));
        count=1;
        for k=pc_list{i}
            place_cells(centroids{i}==k)=1;
            SI_cells(centroids{i}==k)=SI{i}(count);
            count=count+1;
        end

        all_cells=mat2cell(all_cells,winSize/bins.*ones(1,bins),winSize/bins.*ones(1,bins));
        place_cells=mat2cell(place_cells,winSize/bins.*ones(1,bins),winSize/bins.*ones(1,bins));
        SI_cells(SI_cells==0)=NaN;
        SI_cells=mat2cell(SI_cells,winSize/bins.*ones(1,bins),winSize/bins.*ones(1,bins));

        all_cells=cellfun(@(x) sum(sum(x)),all_cells);
        place_cells=cellfun(@(x) sum(sum(x)),place_cells);
        SI_map(map_index==i)=cellfun(@(x) mean(mean(x,'omitnan'),'omitnan'),SI_cells);

        map(map_index==i)=place_cells./all_cells;
        
        mimg_map(mimg_index==i)=mimg_cell{i};
    end
end

map(isnan(map))=0;
SI_map(isnan(SI_map))=0;

figure
imagesc(mimg_map)
pbaspect([numWin(2) numWin(1) 1]);
colormap gray;
xticks(0:winSize/fov*winSize/2:winSize*numWin(2));
xticklabels(split(num2str(0:winSize/2:winSize*numWin(2))));
xlabel('ML (\mum)');
[idx,~]=find(index==bregma);
yticks(unique([fliplr((idx-1)*winSize:-winSize/fov*winSize/2:0) (idx-1)*winSize:winSize/fov*winSize/2:winSize*numWin(1)]));
label=fliplr(0:winSize/2:winSize*numWin(1));
yticklabels(split(num2str(label-label(length(fliplr((idx-1)*winSize:-winSize/fov*winSize/2:0))))));
ylabel('AP (\mum)');

figure
imagesc(Smooth(map,[smoothing smoothing]))
pbaspect([numWin(2) numWin(1) 1]);
colormap jet;
c=colorbar; c.Label.String='Fraction place cells';
xticks(0:bins/fov*winSize/2:bins*numWin(2));
xticklabels(split(num2str(0:winSize/2:winSize*numWin(2))));
xlabel('ML (\mum)');
yticks(unique([fliplr((idx-1)*bins:-bins/fov*winSize/2:0) (idx-1)*bins:bins/fov*winSize/2:bins*numWin(1)]));
label=fliplr(0:winSize/2:winSize*numWin(1));
yticklabels(split(num2str(label-label(length(fliplr((idx-1)*winSize:-winSize/fov*winSize/2:0))))));
ylabel('AP (\mum)');

figure
imagesc(Smooth(SI_map,[smoothing smoothing]))
pbaspect([numWin(2) numWin(1) 1]);
colormap jet;
c=colorbar; c.Label.String='Spatial info (bits)';
xticks(0:bins/fov*winSize/2:bins*numWin(2));
xticklabels(split(num2str(0:winSize/2:winSize*numWin(2))));
xlabel('ML (\mum)');
yticks(unique([fliplr((idx-1)*bins:-bins/fov*winSize/2:0) (idx-1)*bins:bins/fov*winSize/2:bins*numWin(1)]));
label=fliplr(0:winSize/2:winSize*numWin(1));
yticklabels(split(num2str(label-label(length(fliplr((idx-1)*winSize:-winSize/fov*winSize/2:0))))));
ylabel('AP (\mum)');