classdef multiplane < handle
    properties (GetAccess = 'public', SetAccess = 'protected')
        planes      % array of LFP objects
        tree        % tree structure used for dirTree
        ROI = struct('overlap',[])  % overlapping/non-overlapping ROIs
        
        ops = struct('target', [], 'lfp_ops', []);
    end
    
    methods
        function obj = multiplane(path)
            if nargin < 1
                path = pwd;
                disp('No directory specified. Defaulting to current directory.');
            end
            
            root = dir(path);
            temp = strsplit(path, {'\\','/'});
            if isunix
                obj.tree.parent = strjoin(temp(1:end-1), '/');
            elseif ispc
                obj.tree.parent = strjoin(temp(1:end-1), '\\');
            else
                error('What kind of weird operating system are you using?');
            end
            obj.tree.name = temp{end};
            obj.planes = LFP.empty(0);
            count = 1;
            for ii = 3:length(root)
                obj.tree.child(ii-2).name = root(ii).name;
                list = dir(fullfile(path, root(ii).name, '*.abf'));
                for jj = 1:length(list)
                    obj.tree.child(ii-2).child(jj).name = list(jj).name;
%                     obj.planes(count) = LFP({ fullfile(path, root(ii).name, list(jj).name) });
%                     obj.planes(count).enregistrer();
                    obj.planes{count} = LFP.empty(0);
                    count = count + 1;
                end
            end
%             obj.set_target();
        end
        
        function LFP(obj, x)
            obj.planes{x} = LFP([ { fullfile(obj.tree.parent, obj.tree.name, obj.tree.child(x)) } obj.lfp_ops]);
        end
        
        function set_lfp_ops(obj, varargin)
            obj.ops.lfp_ops = varargin;
        end
        
        function set_target(obj, varargin) % define the targets of analysis
            % {'date', 'dates'} - dates to process as an array of datetime
            if ~isfield(obj.ops, 'target')
                obj.ops.target.dates = [];
                obj.ops.target.sessions = NaN;
                obj.ops.target.planes = NaN;
            end
            
            count = 1;
            while count <= length(varargin)
                switch lower(varargin{count})
                    case {'date', 'dates'}
                        obj.ops.target.dates = varargin{count + 1};
                    case {'session', 'sessions'}
                        obj.ops.target.sessions = varargin{count + 1};
                    case {'plane', 'planes'}
                        obj.ops.target.planes = varargin{count + 1};
                    otherwise
                        error(['''' varargin{count} ''' is not a valid parameter']);
                end
                count = count + 2;
            end
            
            obj.ops.target.index = false(length(obj.planes), 1);
            for ii = 1:length(obj.planes)
                obj.ops.target.index(ii) = (ismember(obj.planes{ii}.session.date, obj.ops.target.dates) || isempty(obj.ops.target.dates)) && ...
                                           (ismember(str2double(obj.planes(ii).session.id), obj.ops.target.sessions) || isnan(obj.ops.target.sessions)) && ...
                                           (ismember(obj.planes(ii).twop.plane, obj.ops.target.planes) || isnan(obj.ops.target.planes));
            end
            obj.ops.target.index = find(obj.ops.target.index)';
        end
        
        function choose_days(obj)
            dates = {obj.tree.child.name};
            idx = listdlg('liststring', dates);
        end
        
        function perform_analysis(obj)
            for ii = obj.ops.target.index
                obj.planes(ii).perform_analysis;
                obj.planes(ii).enregistrer();
            end
        end
        
        function intersect(obj) %find overlapping/non-overlapping ROIs
            if isempty(obj.ROI.overlap)
                obj.ROI.overlap = cell(length(obj.planes));
            end
            
            count = 1;
            for ii = obj.ops.target.index
                for jj = obj.ops.target.index
                    if isempty(obj.ROI.overlap{ii,jj})
                        obj.ROI.overlap{ii,jj} = obj.register(ii, jj);
                    end
                    disp(['Done ' num2str(count) '/' num2str(length(obj.ops.target.index)^2)])
                    count = count + 1;
                end
            end
        end
        
        function dir(obj) % print the directory tree
            dirTree(obj.tree);
        end
        
        plot(obj, type);
        
        overlap = register(obj, a, b, varargin);
    end
end