function player(obj)

pos=1;

h=figure;
ax1=subplot(2,1,1);
ax2=subplot(2,1,2);

showFrame;
showTrace;

prev=uicontrol(h,'style','pushbutton','string','previous','position',[20 h.Position(4)-20 60 20],'callback',@prev_frame);
next=uicontrol(h,'style','pushbutton','string','next','position',[h.Position(3)-20-60 h.Position(4)-20 60 20],'callback',@next_frame);

seek_box=uicontrol(h,'style','checkbox','string','seek','callback',@seek);


    function showFrame()
        imagesc(ax1,obj.cam.read(pos))
        axis(ax1,'image')
    end

    function showTrace()
        hold(ax2,'off');
        plot(ax2,obj.traces);
        hold(ax2,'on');
        if ~isempty(obj.mvt)
            plot(ax2,find(obj.mvt),obj.traces(obj.mvt),'.');
        end
        plot(ax2,pos,obj.traces(pos)+5,'vr');
    end

    function seek(src, event)
        while(seek_box.Value)
            [x,~]=ginput(1);
            if(x<0)
                break;
            end
            pos=round(x);
            showFrame;
            showTrace;
        end
    end

    function next_frame(src, event)
        if pos < obj.num_frames
            pos=pos+50;
            showFrame;
            showTrace;
        end
    end
    function prev_frame(src, event)
        if pos > 1
            pos=pos-50;
            showFrame;
            showTrace;
        end
    end
end