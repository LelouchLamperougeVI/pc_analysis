classdef ensemble < LFP

    properties (GetAccess = 'public', SetAccess = 'public')
        order
        pc_order
    end
    properties (GetAccess = 'public', SetAccess = 'private')
        ops
        SCE
        MUA
        clust_SCE
        clust_MUA
        R
        null_R
        tree % agglomerative tree
        h_thres % threshold from clustering
        clust % ensemble clusters
        clust_order % order neurons by ensembles
        swr_stack
        swr_all
        swr_null_all
        swr_t
        swr_class
        colours
    end
    
    methods
        function obj = ensemble(varargin)
            obj@LFP(varargin);
            if isempty(obj.analysis)
                warning('no analysis available');
            else
                obj.pc_order = obj.analysis.order;
            end
        end
        
        function perform_analysis(obj)
            perform_analysis@LFP(obj);
        end
        
        function duration(obj,thres) % remove SCEs outside duration threshold (thres = [lower_lim upper_lim])
            idx=obj.SCE.dur<thres(1) | obj.SCE.dur>thres(2);
            obj.SCE.dur(idx)=[];
            obj.SCE.on(idx)=[];
            obj.SCE.peak(idx)=[];
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
        plot(obj,type,varargin);
    end
    
    methods (Access=private)
        [clust, s, ticks] = silhouette_cluster(obj,Z,D,e_size);
    end
    
end