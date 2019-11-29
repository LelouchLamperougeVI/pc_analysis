classdef dirTree < handle
    properties (GetAccess = 'public', SetAccess = 'protected')
        name % name of current folder
        parent % parent path
        children = dirTree.empty(0) % list of children also of class dirTree
        filters = struct('isdir',NaN, 'regex',{{NaN}})
        fcn = {NaN} % list of functions to perform on the name of the directory/file
        var % parsed variable
    end
    
    methods
        function obj = dirTree(varargin)
            target = obj.parse_inputs(varargin);
            
            [obj.parent, name, ext] = fileparts(target);
            obj.name = strcat(name, ext);
            
            obj.filter;
            if isempty(obj)
                return
            end
            
            subs = dir(target);
            if length(subs)==1 && ~subs(1).isdir
                return
            end
            subs(1:2) = [];
            
            count = 1;
            for ii = 1:length(subs)
                child = dirTree(fullfile(subs(ii).folder, subs(ii).name), 'isdir', obj.filters.isdir(2:end), 're', obj.filters.regex(2:end), 'fcn', obj.fcn(2:end));
                if ~isempty(child)
                    obj.children(count) = child;
                    count = count + 1;
                end
            end
        end
        
        function [nodes, txtlist] = tree(obj, parent, child)
            if nargin<2
                parent = 0;
                child = 1;
            end
            nodes = parent;
            txtlist = {obj.name};
            parent = child;
            for ii = 1:length(obj.children)
                child = length(nodes) + parent;
                [node, txt] = obj.children(ii).tree(parent, child);
                nodes = [nodes node];
                txtlist = [txtlist txt];
            end
            if nargin<2
                [x, y] = treelayout(nodes);
                figure;
                treeplot(nodes);
                for ii = 1:length(txtlist)
                    text(x(ii), y(ii) - .01, txtlist{ii}, 'interpreter','none', 'fontweight','bold', 'backgroundcolor','w', 'edgecolor','k');
                end
                camroll(90);
            end
        end
    end
    
    methods (Access = private)
        function target = parse_inputs(obj, inputs)
            count = 1;
            if mod(length(inputs), 2)
                target = inputs{count};
                count = count + 1;
            else
                target = pwd;
            end
            
            while count <= length(inputs)
                switch lower(inputs{count})
                    case 'isdir'
                        obj.filters.isdir = inputs{count + 1};
                    case {'re', 'regex', 'regexp'}
                        obj.filters.regex = inputs{count + 1};
                    case 'fcn'
                        obj.fcn = inputs{count + 1};
                    otherwise
                        error(['''' inputs{count} ''' is not a valid parameter']);
                end
                count = count + 2;
            end
        end
        
        function filter(obj)
            if isempty(obj.filters.isdir)
                isd = false;
            elseif isnan(obj.filters.isdir(1))
                isd = true;
            else
                isd = isfolder(fullfile(obj.parent, obj.name)) == obj.filters.isdir(1);
            end
            if isempty(obj.filters.regex)
                re = false;
            elseif isnan(obj.filters.regex{1})
                re = true;
            else
                re = regexp(obj.name, obj.filters.regex{1}, 'once');
            end
            if ~( isd && ~isempty(re) )
                obj = obj.empty(0);
                return
            end
            if isa(obj.fcn{1}, 'function_handle')
                obj.var = obj.fcn{1}(obj.name);
            end
        end
    end
end