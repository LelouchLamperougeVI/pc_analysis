function mask = gui_get_bounds(obj)
h=figure;
imagesc(obj.cam.read(1));
hold on
count=1;
uicontrol('style','pushbutton','string','Skip 10 frames','callback',@next_frame);
coor=[];
x=-1;y=-1;
while ~isempty(x)
    if x>0 && y>0
        coor=[coor; [x y]];
        plot(x,y,'rd');
    end
    if size(coor,1)>1
        plot([coor(end-1,1) coor(end,1)], [coor(end-1,2) coor(end,2)], 'b');
    end
    [x,y]=ginput(1);
end

coor=round(coor);
coor=[coor;coor(1,:)];
mask=false(obj.dims(1),obj.dims(2));
for i=1:size(coor,1)-1
    [seg,idx]=max(abs(diff(coor(i:i+1,:))));
    temp=round(linspace(coor(i,mod(idx,2)+1), coor(i+1,mod(idx,2)+1), seg+1));
    seg=zeros(seg+1,2);
    seg(:,idx)=linspace(coor(i,idx), coor(i+1,idx), length(temp));
    seg(:,mod(idx,2)+1)=temp;
    mask(sub2ind(size(mask),seg(:,1),seg(:,2)))=true;
end

mask=imfill(mask,'holes');

[x,y]=find(mask);
plot(x,y,'*c');

    function next_frame(src, event)
        count=count+10;
        if count<=obj.num_frames
            hold off
            imagesc(obj.cam.read(count));
            hold on
        end
    end
end
