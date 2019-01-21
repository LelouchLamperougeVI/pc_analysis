classdef basic_ensemble < handle
% A very basic analysis for detecting synchronous calcium events (SCE) by
% performing a sum of Q over zscored deconv
% This should serve as a temporary solution for the Philosophical
% Transactions paper, while we work on finalizing the MI method
% Construction:
%   deconv
%   ts:             timestamps for 2p frames in seconds
%
%   Name, Value:
%   'sig':          gaussian smoothing kernel SD (default - 50 ms)
%   'thres':        threshold for SCE expressed in SDs (default - 3 SD)
%   'off_thres':    threshold for SCE onset/offset (default - 1 SD)
%   'gaps':         distance in seconds for 2 SCEs to be considered as a singular event (default - 250 ms)
%
% Output: 
%   SCE.
%           dur:       duration of a SCE event
%           on:        onset time of SCEs
%   MUA:    means multi-unit dF/F

    properties (GetAccess = 'public', SetAccess = 'public')
        order
        lfp % object of type LFP
    end
    properties (GetAccess = 'public', SetAccess = 'private')
        ops
        deconv
        ts
        SCE
        MUA
        spec = struct('spectrum',[],'f',[],'t',[],'norm',[])
        R
        null_R
        tree % agglomerative tree
        clust % ensemble clusters
    end
    
    methods
        function obj = basic_ensemble(deconv,ts,varargin)
            obj.deconv=deconv;
            obj.ts=ts;
            obj.detect_sce(varargin);
        end
        
        function duration(obj,thres) % remove SCEs outside duration threshold (thres = [lower_lim upper_lim])
            idx=obj.SCE.dur<thres(1) | obj.SCE.dur>thres(2);
            obj.SCE.dur(idx)=[];
            obj.SCE.on(idx)=[];
            obj.SCE.peak(idx)=[];
        end
        
        set_ops(obj,varargin);
        
        detect_sce(obj,varargin);
        sce_spectrum(obj,varargin);
        corr(obj);
        hclust(obj);
        plot(obj,type);
    end
    
%     methods (Access = private)
%         set_ops(obj,inputs);
%     end
end