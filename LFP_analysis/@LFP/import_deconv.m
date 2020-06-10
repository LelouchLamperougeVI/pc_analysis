function import_deconv(obj,deconv) % import deconv
if any(deconv(:) < 0)
    warning('Deconv contains negative firing rates. These would be set to 0.');
    deconv(deconv < 0) = 0;
end

if exist(fullfile(obj.session.wd, 'fkd18k'), 'file')
    fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'r');
    if fscanf(fid, '%d')
        obj.twop.deconv = stupid_windows_fs(deconv);
        disp('18k correction applied');
    else
        obj.twop.deconv = deconv;
        if size(deconv,1)>18000
            disp('18k correction unnecessary');
        end
    end
    fclose(fid);
    return
end
if size(deconv,1)>18000 && ~isempty(obj.behavior)
    tcs.tt=obj.twop.ts';
    [vel,temp] = convert_behavior(obj.behavior,tcs,deconv);
    vel = vel.unit_vel;
    temp = fast_smooth(temp,obj.twop.fs*.5);
    if size(temp,1) < length(vel); vel=vel(1:size(temp,1)); end
    ori_r = corr(temp, vel');
    
    [vel,temp] = convert_behavior(obj.behavior,tcs,stupid_windows_fs(deconv));
    vel = vel.unit_vel;
    temp = fast_smooth(temp,obj.twop.fs*.5);
    if size(temp,1) < length(vel); vel=vel(1:size(temp,1)); end
    alt_r = corr(temp, vel');
    
    if prctile(abs(alt_r), 70) > prctile(abs(ori_r), 70)
        disp('DANGER! A rupture in Space-Time continuum has been detected!');
        uinp = input('Should I fix it for you? Y/N [Y] ', 's');
        if strcmpi(uinp,'y') || strcmpi(uinp,'yes') || isempty(uinp)
            obj.twop.deconv=stupid_windows_fs(deconv);
            disp('OK. I''ve fixed it for you. But know that this only works for behavioural data. It''s up to you to remember to correct it during rest.');
            fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'w');
            fprintf(fid, '%d', 1);
            fclose(fid);
            return
        end
        disp('I surely hope you know what you''re doing...');
    end
end
obj.twop.deconv=deconv;
fid = fopen(fullfile(obj.session.wd, 'fkd18k'), 'w');
fprintf(fid, '%d', 0);
fclose(fid);