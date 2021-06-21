clear all

traj_thres = .2;
root = '/mnt/md0/Data/RSC_M2';
animals = dir(fullfile(root, 'RSC*'));
animals = {animals.name};

clust_stacks = cell(2,1);
trajectories = cell(2,1);
loc_clust = cell(2,1);
clusts = cell(2,1);
pev = cell(2, 1);
for a = 1:length(animals)
    sessions = dir(fullfile(root, animals{a}));
    sessions = {sessions.name};
    sessions = sessions( ~cellfun(@isempty, regexp(sessions, '^\d\d\d\d_\d\d_\d\d')) );
    
    for s = 1:length(sessions)
        rest1 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_1.abf']));
        rest1.set_ops('e_size',5);
        rest1.set_ops('clust_method','thres');
        rest1.set_ops('sig', .2);
        rest1.remove_mvt;
        rest1.cluster;
        
        rest2 = ensemble(fullfile(root, animals{a}, sessions{s}, [sessions{s} '_3.abf']));
        rest2.set_ops('e_size',5);
        rest2.set_ops('clust_method','thres');
        rest2.set_ops('sig', .2);
        rest2.remove_mvt;
        rest2.cluster;
        
        % v checkpoint
    
        clusts{1} = cat(2, clusts{1}, {rest1.ensembles.clust});
        clusts{2} = cat(2, clusts{2}, {rest2.ensembles.clust});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        
        analysis = rest1.analysis;
        L = length(analysis.Pi);
        thres=noRun(analysis.behavior.unit_vel);
        thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
        n = analysis.original_deconv;
        dt = .3; % 300 milliseconds
        dt = round(analysis.fs * dt);
        x = analysis.behavior.unit_pos;
            
        for c = 1:length(rest1.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest1.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            clust_stacks{1} = cat(1, clust_stacks{1}, {stack});
            trajectories{1} = cat(2, trajectories{1}, any(stack > traj_thres, 2));
            
            try
                md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false);
                pev{1} = cat(1, pev{1}, {[md.pc.pev(:), md.cc.pev(:)]});
            catch
                pev{1} = cat(1, pev{1}, {[]});
            end
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{1} = cat(1, loc_clust{1}, {temp});
        
        g = cellfun(@(x) zeros(size(x, 1), 1), rest1.analysis.width, 'uniformoutput', false);
        
        analysis = rest2.analysis;
        L = length(analysis.Pi);
        thres=noRun(analysis.behavior.unit_vel);
        thres=(analysis.behavior.unit_vel>thres | analysis.behavior.unit_vel<-thres) & (analysis.behavior.trials(1) < analysis.behavior.frame_ts & analysis.behavior.trials(end) > analysis.behavior.frame_ts);
        n = analysis.original_deconv;
        dt = .3; % 300 milliseconds
        dt = round(analysis.fs * dt);
        x = analysis.behavior.unit_pos;
        
        for c = 1:length(rest2.ensembles.clust)
            list = intersect(rest1.analysis.pc_list, rest2.ensembles.clust{c});
            stack = rest1.analysis.stack(:, list);
            clust_stacks{2} = cat(1, clust_stacks{2}, {stack});
            trajectories{2} = cat(2, trajectories{2}, any(stack > traj_thres, 2));
            
            try
                md = pc_cc_simanneal(n(:, list), x, dt, 'reject', ~thres, 'prog', false);
                pev{2} = cat(1, pev{2}, {[md.pc.pev(:), md.cc.pev(:)]});
            catch
                pev{2} = cat(1, pev{2}, {[]});
            end
            
            for l = 1:length(list)
                g{list(l)} = c .* ones(size(rest1.analysis.width{list(l)}, 1), 1);
            end
        end
        temp = cell2mat({g{~cellfun(@isempty, rest1.analysis.width)}}');
        loc_clust{2} = cat(1, loc_clust{2}, {temp});
    end
end
