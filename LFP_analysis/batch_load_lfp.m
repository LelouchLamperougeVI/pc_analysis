function [lfp,frame_ts,fs]=batch_load_lfp(date,sessions,ch)

if nargin<3
    ch=6;
end

lfp=cell(1,length(sessions));
frame_ts=lfp;
count=1;
for i=sessions
    [signal,fs]=load_lfp([ch 1],['\' date '_' num2str(i) '.abf'], pwd);
    lfp{count}=signal(:,1);
    signal=signal(:,2);
    signal=signal<3;
    signal=find(get_head(signal));
    if signal(1)==1; signal(1)=[]; end
    frame_ts{count}=signal./fs;
    count=count+1;
end