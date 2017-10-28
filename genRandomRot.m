%% Generate random rotation trials for shuffle test, because stupid me used constant velocity -.-

function [unit_rot,trials] = genRandomRot(ts)
    maxV=60;
    minV=10;
    dT=diff(ts);
    dT(end+1)=dT(end);
    unit_rot=zeros(1,length(ts));
    pre=1;
    trials=[];
    while(1)
        direction=ceil(2*rand(1))-1;
        if direction==0
            direction=-1;
        end
        vel=direction*((maxV-minV)*(rand(1))+minV);
        
        try
            next=find(ts>=(ts(pre)+360/abs(vel)),1);
        catch
            return
        end
        if isempty(next)
            unit_rot(pre:end)=dT(pre:end).*vel;
            unit_rot(pre:end)=cumsum(unit_rot(pre:end));
            if vel>0
                unit_rot(pre:end)=unit_rot(pre:end)-unit_rot(pre);
            else
                unit_rot(pre:end)=unit_rot(pre:end)-unit_rot(end);
            end
            return
        end
        trials=[trials ts(next)];
        unit_rot(pre:next)=dT(pre:next).*vel;
        unit_rot(pre:next)=cumsum(unit_rot(pre:next));
        if vel>0
            unit_rot(pre:next)=unit_rot(pre:next)-unit_rot(pre);
        else
            unit_rot(pre:next)=unit_rot(pre:next)-unit_rot(next);
        end
        pre=next+1;
    end