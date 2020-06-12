%% A temporary fix for compatibility with tablet
% Still in the process of figuring out how to process behavioral data from
% tablet. In the mean time, this can be used to extract raw and deconv signals.
% N.B. it saves automatically you big dummy

type='Place cells';
type='survey';
% type='HD cells';
newdir={'rsc034\2017_07_14\1','rsc032\2017_07_14\1','rsca1002\2017_07_14\1','rsca1002\2017_07_14\2'};
newdir={'ca1012\2017_07_10\1','ca1013\2017_07_14\1'};
newdir={'rsc032\2017_07_26\1','rsc032\2017_07_26\2'};
% newdir={'rsc034\2017_07_25\1','rsc034\2017_07_25\2','rsc034\2017_07_25\3'};
% newdir={'ca1012\2017_07_26\1','ca1012\2017_07_26\2'};
newdir={'rsc034\2017_07_26\1','rsc034\2017_07_26\2'};
newdir={'rsc032\2017_07_29\1','rsc032\2017_07_29\2','rsc032\2017_07_29\3'};
newdir={'rsc032\2017_07_21\1'};
newdir={'ca1013\2017_08_11\2'};
newdir={'rsc032\2017_07_22\1','rsc032\2017_07_22\2','rsc032\2017_07_22\3','rsc032\2017_07_22\4'};
newdir={'rsc034\2017_09_13\1','rsc034\2017_09_13\2','rsc034\2017_09_13\3','rsc034\2017_09_13\4','rsc034\2017_09_13\5'};
newdir={'ca1011\2017_10_04\2'};
newdir={'ca1011\2017_10_07\2'};
parfor i=1:length(newdir)
    expt = frGetExpt_leth(newdir{i});
    expt = ops_update_xml(expt);

    subdirs = dir(expt.dirs.reggreenpn);
    okdir = [subdirs(:).isdir];
    namesubdirs = {subdirs(okdir).name}';
    namesubdirs(ismember(namesubdirs, {'.', '..'})) = [];

    newdirec = fullfile(newdir{i}, namesubdirs{1});
    expt = frGetExpt_leth2(newdirec);
    expt = ops_update_xml(expt);
    MakeDirCD(expt);
%     MakeDirCD(newdir{i});

    GetMaskSuite2P(expt);

    tcs_stack(i) = frExtractTimeCourses(expt,'Overwrite',true,'UpdateRaw',true,'AdaptBaseline',false,'RingSubtract',true,'MaskOp','manual','MaskAlign','aligned');
%     tcs = frExtractTimeCourses(expt,'Overwrite',true,'UpdateRaw',true,'AdaptBaseline',false,'RingSubtract',true,'MaskOp','auto','MaskAlign','aligned');

    [c,~,~,~,~,sp] = constrained_foopsi(double(tcs(i).ratio));
%     ratio_model = reshape(c, size(tcs(i).ratio,1), size(tcs(i).ratio,2));
    deconv_cell{i} = reshape(sp, size(tcs(i).ratio,1), size(tcs(i).ratio,2));
end

% Uncomment if you're too lazy to segment by hand (test phase)
% for i=1:length(newdir)
%     expt = frGetExpt_leth(newdir{i});
%     expt = ops_update_xml(expt);
%     expt = frGetExpt_leth2(newdirec);
%     expt = ops_update_xml(expt);
%     load([expt.dirs.analrootpn '\masks_neurons.mat']);
%     tt(i)=cell_filt(tt(i),maskNeurons);
% end

for i=1:length(newdir)
    cd(['Y:\homes\haoran.chang\data\' type]);
    split=strsplit(newdir{i},'\');
    try
        mkdir(split{2});
    catch
    end
    cd(split{2});
    
    tcs=tcs_stack(i);
    deconv=deconv_cell(i);
    
    save([strjoin(strsplit(newdir{i},'\'),'_') '.mat'], 'tcs', 'deconv');
end
