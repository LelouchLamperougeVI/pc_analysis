function obj = set(obj, ops)
% Set options for LFP object
% If ops = []; reset properties to default

% Channels = [1 2 3 5 5 5]'; %old behavior for RSC RRR (phil trans B)
% Channels = [1 2 3 6 nan nan]'; %old behavior for RSC RRR (phil trans B)
% Channels = [1 2 3 6 nan nan 5]'; % licks
% Channels = [1 2 3 5 nan nan nan]'; % Dun topography collab
% Channels = [1 2 3 5 nan nan]'; %old behavior for RSC RRR (phil trans B)
% Channels = [1 2 3 5 5 6]'; %old old behavior for RSC RRR, also works for Ingrid's RRR
% Channels = [1 2 3 6 8 7]'; %new behavior for RSC RRR
% Channels = [1 2 3 4 8 7]'; %new vr behavior for RSC RRR
% Channels = [1 3 2 6 5 5]'; % lesion across days
% Channels = [1 2 3 4 5 5]'; % aubrey's data
% Channels = [1 4 5 6 nan nan]'; % aubrey's data
% Channels = [1 2 3 5 5 5]'; % haoran's data

if isempty(ops) % reset default properties
    obj.abf.Channels = [1 2 3 5 nan nan nan];
    obj.abf.ChannelNames = {'2p';'chA';'chB';'rwd';'cam';'lfp';'lck'};
    
    obj.lfp.f60_env = 1;
    obj.lfp.ops.nfft = 2^16;
    obj.lfp.ops.down_fs = 2e3;
    obj.lfp.ops.spec_wdw = .5;
    obj.lfp.ops.swr_cyc = 3;
    obj.lfp.ops.swr_gap = .25;
    
    obj.topo.FOV = 835.76 .* ones(2,1);
    
    obj.twop.planes.planes = 1;
    obj.twop.planes.numplanes = 1;
    obj.twop.planes.plane_names = {'single plane'};
    obj.twop.planes.ol = .8;
    obj.twop.planes.maxStep = 100;
    obj.twop.planes.maxR = .7;
    
    obj.intern.fn = [];
    obj.intern.path = [];
    obj.intern.op_mode = 'abf';
end

count = 1;
while count <= length(ops)
    switch lower(ops{count})
        case {'ch', 'chan', 'channel', 'channels'}
            obj.abf.Channels = ops{count + 1};
            obj.update_channels;
        case {'plane', 'planes'}
            obj.twop.planes.planes = ops{count+1};
        case {'ol', 'overlap'}
            obj.twop.planes.ol = ops{count+1};
        case 'maxstep'
            obj.twop.planes.maxStep = ops{count+1};
        case 'maxr'
            obj.twop.planes.maxR = ops{count+1};
        otherwise
            exp = regexp(ops{count}, {'\.abf$', 'behavior\.mat$'}, 'once');
            idx = ~cellfun(@isempty, exp);
            if any(idx)
                obj.intern.path = strsplit(ops{count}, {'\\','/'});
                obj.intern.fn = obj.intern.path{end};
                if isunix
                    obj.intern.path = strjoin( obj.intern.path(1:end-1), '/');
                elseif ispc
                    obj.intern.path = strjoin( obj.intern.path(1:end-1), '\');
                end
                count = count - 1;

                switch find(idx)
                    case 1
                        obj.intern.op_mode = 'abf';
                    case 2
                        obj.intern.op_mode = 'beh';
                end
            else
                error(['The argument ''' ops{count} ''' is not a valid property for the class ' class(obj)]);
            end
    end
    
    count = count + 2;
end