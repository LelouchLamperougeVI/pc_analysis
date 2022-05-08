classdef ensemble < LFP
    properties (GetAccess = 'public', SetAccess = 'private')
        ops
        ensembles = struct('SCE', [], 'MUA', [], 'R', [], 'null_R', [], 'clust', [])
        hiepi
        bayes
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        intern
    end
    
    methods
        function obj = ensemble(varargin)
            obj@LFP(varargin);
            if isempty(obj.analysis)
                warning('no analysis available');
            end
        end
        
        function duration(obj,thres) % remove SCEs outside duration threshold (thres = [lower_lim upper_lim])
            idx=obj.ensembles.SCE.dur<thres(1) | obj.ensembles.SCE.dur>thres(2);
            obj.ensembles.SCE.dur(idx)=[];
            obj.ensembles.SCE.on(idx)=[];
            obj.ensembles.SCE.peak(idx)=[];
        end
        
        function cluster(obj)
%             obj.duration(obj.ops.sce_dur);
            obj.hclust;
        end
        
        set_ops(obj,varargin);
        
        detect_sce(obj,varargin);
        sce_spectrum(obj,varargin);
        corr(obj);
        hclust(obj);
        [stack,all,t]=swr_window(obj);
        topography(obj);
        make_colours(obj);
        classi_swr(obj, wdw, p);
        [decoded, P, pos, err] = bayes_infer(obj, varargin);
        [decoded, P, pos, err] = bayes_infer2(obj, varargin);
        plot(obj,type,varargin);
        hiepi = hPICA(obj);
    end
    
    methods (Access=private)
        [clust, s, ticks] = silhouette_cluster(obj,Z,D,e_size);
    end
    
end