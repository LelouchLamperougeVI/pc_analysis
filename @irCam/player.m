function player(obj)

pos=1;

h=figure;
ax1=subplot(2,1,1);
ax2=subplot(2,1,2);
xlabel(ax2,'frame');
ylabel(ax2,'pixel change');

showFrame;
showTrace;
fullscale;

prev=uicontrol(h,'style','pushbutton','string','previous','position',[20 h.Position(4)-20 60 20],'callback',@prev_frame);
next=uicontrol(h,'style','pushbutton','string','next','position',[h.Position(3)-20-60 h.Position(4)-20 60 20],'callback',@next_frame);

full=uicontrol(h, 'style','pushbutton', 'string','full scale', 'position',[h.Position(3)/2 20 60 20], 'callback',@fullscale);

seek_box=uicontrol(h,'style','checkbox','string','seek','callback',@seek);


    function showFrame()
        imagesc(ax1,obj.cam.read(pos))
        axis(ax1,'image')
    end

    function showTrace()
        px = get(ax2,'XLim');
        py = get(ax2,'YLim');
        
        hold(ax2,'off');
        plot(ax2,abs(obj.traces));
        hold(ax2,'on');
        if ~isempty(obj.mvt)
            plot(ax2,find(obj.mvt==1), -1.*ones(sum(obj.mvt==1),1),'r.', 'markersize',5);
            plot(ax2,find(obj.mvt==2), -1.*ones(sum(obj.mvt==2),1),'g.', 'markersize',5);
        end
        plot(ax2,[pos pos],[-10 10],'r');
        
        set(ax2, 'XLim',px);
        set(ax2, 'YLim',py);
    end

    function seek(src, event)
        while(seek_box.Value)
            [x,~,mousebutt]=ginput(1);
            if(x<0)
                break;
            end
            if mousebutt == 1
                pos=round(x);
                showFrame;
                showTrace;
            elseif mousebutt == 3
                toggle_mvt(round(x));
                showTrace;
            end
        end
    end

    function toggle_mvt(x)
        if ~obj.mvt(x)
            return
        end
        heads = find( get_head(logical(obj.mvt)) );
        tails = get_head(logical(obj.mvt(end:-1:1)));
        tails = find( tails(end:-1:1) );
        
        [mh,head] = min(abs(x - heads));
        [mt,tail] = min(abs(tails - x));
        
        [~,idx] = min( [mh mt] );
        if idx==1
            idx=head;
        else
            idx=tail;
        end
        idx = heads(idx):tails(idx);
        if all(obj.mvt(idx)==1)
            obj.mvt(idx)=2;
        elseif all(obj.mvt(idx)==2)
            obj.mvt(idx)=1;
        else
            error('whoops, something broke while trying to toggle movement status');
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

    function fullscale(src, event)
        set(ax2, 'XLim',[0 size(obj.traces,1)]);
        set(ax2, 'YLim', [min(obj.traces(:)) max(obj.traces(:))]);
    end
end