% Batch analysis for reactivation
% Context:      Phil. Trans. B 2019 manuscript prelim analysis
% Data folder:  Adam's old computer (AKA my computer) E:\HaoRan\RRR\ {RSC036-38}

% list = {'RSC036', 'RSC037', 'RSC038'};
% list = {'EE001_new', 'PCH017'};
% list = {'PCH017'};
% list = {'RSC032'};
list = {'RSC037'};

%% Batch analysis

for f = 1:length(list)
    
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        clear lfp;
        full = fullfile(root(i).folder, root(i).name);
        
        load(fullfile(full, '2', 'Plane1', 'deconv.mat'));
        lfp = LFP(fullfile(full, [root(i).name '_2.abf']));
        lfp.import_deconv(deconv);
        lfp.perform_analysis;
        analysis = lfp.analysis;
        save(fullfile(full, 'analysis.mat'), 'analysis');
        save(fullfile(full, 'lfp2.mat'), 'lfp');
        disp(['Got ' num2str(length(analysis.pc_list)) ' place cells out of ' num2str(size(deconv,2))])
        
        clear lfp;
        clear ass
        load(fullfile(full, '1', 'Plane1', 'deconv.mat'));
        lfp = LFP(fullfile(full, [root(i).name '_1.abf']));
        lfp.import_deconv(deconv);
        lfp.import_analysis(analysis);
        lfp.remove_mvt;
        lfp.detect_sce;
        ass = lfp.ensemble;
        save(fullfile(full, 'lfp1.mat'), 'lfp');
        save(fullfile(full, 'ass1.mat'), 'ass');
        
        clear lfp;
        clear ass
        load(fullfile(full, '3', 'Plane1', 'deconv.mat'));
        lfp = LFP(fullfile(full, [root(i).name '_3.abf']));
        lfp.import_deconv(deconv);
        lfp.import_analysis(analysis);
        lfp.remove_mvt;
        lfp.detect_sce;
        ass = lfp.ensemble;
        save(fullfile(full, 'lfp3.mat'), 'lfp');
        save(fullfile(full, 'ass3.mat'), 'ass');
    end
    
end


%% Plot correlation

for f = 1:length(list)
    root = dir(list{f});
    root(1:2) = [];
    
    for i = 1:length(root)
        clear lfp;
        full = fullfile(root(i).folder, root(i).name);
        
        if exist(fullfile(full, 'ass3.mat'),'file')
            load(fullfile(full, 'lfp2.mat'));
            order=get_order(lfp.analysis);
            
            load(fullfile(full, 'lfp1.mat'));
            figure;
            subplot(1,3,1);
            deconv=lfp.deconv(:,order);
            deconv(any(isnan(deconv),2),:)=[];
            deconv=ca_filt(deconv);
            deconv=zscore(deconv);
            deconv=fast_smooth(deconv,2);
            imagesc(corr(deconv));
            axis square;
            colormap jet
            
            load(fullfile(full, 'lfp2.mat'));
            subplot(1,3,2);
            deconv=lfp.deconv(:,order);
            deconv(any(isnan(deconv),2),:)=[];
            deconv=ca_filt(deconv);
            deconv=zscore(deconv);
            deconv=fast_smooth(deconv,2);
            imagesc(corr(deconv));
            axis square;
            colormap jet
            
            load(fullfile(full, 'lfp3.mat'));
            subplot(1,3,3);
            deconv=lfp.deconv(:,order);
            deconv(any(isnan(deconv),2),:)=[];
            deconv=ca_filt(deconv);
            deconv=zscore(deconv);
            deconv=fast_smooth(deconv,2);
            imagesc(corr(deconv));
            axis square;
            colormap jet
            
            title(['Mouse :' list{f} ' Date: ' root(i).name]);
        end

    end
    
end