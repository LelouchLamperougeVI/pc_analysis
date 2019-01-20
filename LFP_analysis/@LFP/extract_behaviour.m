function extract_behaviour(obj)

signal=obj.raw;
nSmooth = 1 / (obj.si * 1e-6);

frame_ts=signal(:,1)<1;
frame_ts=get_head(frame_ts);
frame_ts=find(frame_ts);
frame_ts(1)=[];

rwd=signal(:,obj.get_channel('rwd'));
rwd=find(get_head(rwd<1))'; %returns a single index of value 1 for rest trials

ch_a=signal(:,obj.get_channel('chA'));
ch_b=signal(:,obj.get_channel('chB'))';
ch_a=get_head(ch_a>1)';

forward=ch_a&(ch_b<1);
backward=ch_a&(ch_b>1);
dist=forward+backward.*-1;
pos_raw=dist*pi/5*.1;

speed_raw=fast_smooth((pos_raw./median(diff(frame_ts)))',nSmooth);
speed_raw=speed_raw(frame_ts);

if length(rwd)<3
    obj.behavior.speed_raw=speed_raw;
    return
end

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
ts(ts<rwd(1)|ts>rwd(end))=[];
ts=sort([ts rwd]);
pos_raw=pos_raw(ts);
ts=(ts-1).*(obj.si*1e-6);
idx=find(get_head(pos_raw'<5));
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
behavior.speed_raw=speed_raw;
behavior.pos_cum = pos_cum;
behavior.pos_raw = pos_raw;

obj.behavior=behavior;