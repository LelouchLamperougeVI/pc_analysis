classdef RRR < handle
    
    
    properties (GetAccess = 'public', SetAccess = 'protected')
        lfp = []
        EV
    end
    
    methods
        function obj = RRR(path)
            if nargin < 1
                path = uigetdir(pwd, 'select LFP objects root dir');
            end
            
            list = dir(path);
            idx = arrayfun(@(x) regexpi(x.name, 'lfp\d.mat'), list, 'uniformoutput',false);
            idx = cellfun(@isempty, idx);
            list = {list(~idx).name};
            if length(list)~=3
                error('You must have exactly 3 lfp objects');
            end
            
            for ii = 1:length(list)
                lfp = load(fullfile(path, list{ii}));
                try lfp = lfp.lfp; catch error(['non-existent field ''lfp'' in ' list{ii}]); end
%                if ~isa(lfp,'LFP') && ~isa(lfp,'ensemble'); error([list{ii} ' does not contain a valid LFP class or subclass']); end
                obj.lfp = [obj.lfp lfp];
            end
        end
        
        function [EV, REV] = explained_variance(obj)
            [EV, REV] = ev(obj.lfp(1), obj.lfp(2), obj.lfp(3));
            obj.EV.ev = EV;
            obj.EV.rev = REV;
        end
    end
end