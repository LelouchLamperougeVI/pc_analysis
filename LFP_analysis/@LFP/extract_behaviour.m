function extract_behaviour(obj)

signal=obj.raw;
nSmooth = 1 / (obj.si * 1e-6);

frame_ts=signal(:,obj.get_channel('2p'))<1;
frame_ts=get_head(frame_ts);
frame_ts=find(frame_ts);
frame_ts(1)=[];

if ~isnan(obj.get_channel('lck'))
    licks = signal(:, obj.get_channel('lck')) > 1;
    licks = get_head(licks);
    licks = find(licks);
    licks = knnsearch(frame_ts, licks);
end

rwd=signal(:,obj.get_channel('rwd'));
if median(rwd) > 2.5
    rwd=find(get_head(rwd<1))';
else
    rwd=find(get_head(rwd>1))';
end

ch_a=signal(:,obj.get_channel('chA'));
ch_b=signal(:,obj.get_channel('chB'))';
ch_a=get_head(ch_a>1)';

forward=ch_a&(ch_b<1);
backward=ch_a&(ch_b>1);
dist=forward+backward.*-1;
pos_raw=dist*pi/5*.1;

speed_raw=fast_smooth((pos_raw./median(diff(frame_ts)))',nSmooth);
speed_raw=speed_raw(frame_ts);

behavior.speed_raw_noSmooth = cumsum(pos_raw);
behavior.speed_raw_noSmooth = [0 diff(behavior.speed_raw_noSmooth(frame_ts))./diff(frame_ts)'];
behavior.speed_raw=speed_raw;
if length(rwd)<5
    ts=sort([find(dist~=0) 1 length(dist)]);
    rwd=[];
%     obj.behavior=behavior;
%     return
end

if isempty(rwd)
    pos_raw = cumsum(pos_raw);
else
    count=1;
    while(count<length(rwd))
        if sum(pos_raw(rwd(count):rwd(count+1)))<20
            rwd(count+1)=[];
        else
            count=count+1;
        end
    end
    pos_raw(1:rwd(1)-1)=cumsum(pos_raw(1:rwd(1)-1));
    for i=1:length(rwd)-1
        pos_raw(rwd(i):rwd(i+1)-1)=cumsum(pos_raw(rwd(i):rwd(i+1)-1));
    end
    pos_raw(rwd(end):end)=cumsum(pos_raw(rwd(end):end));

    ts=find(dist~=0);
%     ts(ts<rwd(1)|ts>rwd(end))=[];
    ts=sort([ts rwd 1 length(dist)]);

    spurious = arrayfun(@(x) max(pos_raw(rwd(x):rwd(x+1))), 1:length(rwd)-1); % to detect spots where the valve didn't open or got stuck...
    spur_range = median(spurious);
    spurious = find(spurious > ( spur_range + 20));
    if range( pos_raw(rwd(end) : end) ) > ( spur_range + 20)
        spurious = [spurious length(rwd)];
        rwd = [rwd length(pos_raw)];
    end

    for i = 1:length(spurious)
        temp = pos_raw(rwd(spurious(i)) : rwd(spurious(i) + 1));
        idx = round(range(temp) / spur_range);
        idx = knnsearch(temp', range(temp) .* (1:(idx-1))' ./ idx);

        for j = 1:length(idx)-1
            temp(idx(j):idx(j+1)-1) = temp(idx(j):idx(j+1)-1) - temp(idx(j));
        end
        if diff(temp(end-1:end)) < 0
            temp(idx(end):end-1) = temp(idx(end):end-1) - temp(idx(end));
        else
            temp(idx(end):end) = temp(idx(end):end) - temp(idx(end));
        end

        pos_raw(rwd(spurious(i)) : rwd(spurious(i) + 1)) = temp;
%         rwd = sort( [rwd ( idx' + rwd(spurious(i)) - 1 )] );
    end
end

pos_raw=pos_raw(ts);
ts=(ts-1).*(obj.si*1e-6);
idx=find(get_head(pos_raw'<5));
idx(find(pos_raw(idx(2:end)-1)-pos_raw(idx(2:end))<5)+1)=[];
pos_norm=pos_raw;
pos_cum=pos_raw;
trial=zeros(size(pos_raw));
for i=1:length(idx)-1
    try pos_cum(idx(i):idx(i+1)-1)=pos_cum(idx(i):idx(i+1)-1)+pos_cum(idx(i)-1); catch; end
    pos_norm(idx(i):idx(i+1))=pos_norm(idx(i):idx(i+1))./range(pos_norm(idx(i):idx(i+1)));
    trial(idx(i):idx(i+1))=i;
end
pos_cum(end)=pos_cum(end)+pos_cum(end-1);
pos_norm(idx(end):end)=pos_norm(idx(end):end)./range(pos_norm(idx(end):end));
speed=diff(pos_cum)./diff(ts);
speed=[speed(1) speed];
speed([idx(2:end) idx(2:end)+1])=speed([idx(2:end)-1 idx(2:end)-1]);


behavior.ts = ts;
behavior.pos_norm = pos_norm;
behavior.trial = trial;
behavior.speed = speed;
behavior.pos_cum = pos_cum;
behavior.pos_raw = pos_raw;

obj.behavior=behavior;