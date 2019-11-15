classdef multiplane < handle
    properties (GetAccess = 'public', SetAccess = 'protected')
        planes      % array of LFP objects
        tree        % tree structure used for dirTree
        ROI = struct('overlap',[])  % overlapping/non-overlapping ROIs
        
        ops
    end
    
    methods
        function obj = multiplane(path)
            if nargin < 1
                path = pwd;
                disp('No directory specified. Defaulting to current directory.');
            end
            
            root = dir(path);
            temp = strsplit(path, {'\\','/'});
            obj.tree.name = temp{end};
            obj.planes = LFP.empty(0);
            count = 1;
            for ii = 3:length(root)
                obj.tree.child(ii-2).name = root(ii).name;
                list = dir(fullfile(path, root(ii).name, '*.abf'));
                for jj = 1:length(list)
                    obj.tree.child(ii-2).child(jj).name = list(jj).name;
                    obj.planes(count) = LFP({ fullfile(path, root(ii).name, list(jj).name) });
                    count = count + 1;
                end
            end
        end
        
        function set_target(obj, varargin) % define the targets of analysis
            obj.ops.target.range = []; % range of dates to analyze
            obj.ops.target.sessions = NaN;
            obj.ops.target.planes = NaN;
            
            count = 1;
            while count <= length(varargin)
                switch lower(varargin{count})
                    case 'range'
                        obj.ops.target.range = [varargin{count + 1} varargin{count + 2}];
                        count = count + 1;
                    case 'sessions'
                        obj.ops.target.sessions = varargin{count + 1};
                    case 'planes'
                        obj.ops.target.planes = varargin{count + 1};
                end
                count = count + 1;
            end
        end
        
        function perform_analysis(obj, overwrite)
            if nargin < 2
                overwrite = false;
            end
            for ii = 1:length(obj.planes)
                if isempty(obj.planes(ii).analysis) || overwrite
                    obj.planes(ii).perform_analysis;
                end
            end
        end
        
        function intersect(obj) %find overlapping/non-overlapping ROIs
            obj.ROI.overlap = cell(length(obj.planes));
            
            count = 1;
            for ii = 1:length(obj.planes)
                for jj = 1:length(obj.planes)
                    obj.ROI.overlap{ii,jj} = obj.register(ii, jj);
                    disp(['Done ' num2str(count) '/' num2str(length(obj.planes)^2)])
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