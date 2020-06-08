classdef TOPO < ensemble
    properties (GetAccess = 'public', SetAccess = 'public')
    end
    properties (GetAccess = 'public', SetAccess = 'private')
    end
    
    methods
        function obj = TOPO(varargin)
            obj@ensemble(varargin);
        end
    end
end