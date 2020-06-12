function tcs = cell_filt(tcs,mask)
% Clustering method for identification of noise ROIs
% I know I'm going to forget about this so let me describe the methods
% used:
%
% 1) Compute the autocorrelation for smoothed raw signal of each cell
% 2) Get mean correlation scores between 100 to 400 time lags
% 3) Get peaks and fit to gamma dist
% 4) Cluster features
epsilon=0.05;
MinPts=20;
satisfied=0;
while(~satisfied)
    corr=zeros(tcs.nSamples,tcs.nCells);
    signal=zeros(size(tcs.ratio));
    for i=1:tcs.nCells
        signal(:,i)=smooth(tcs.ratio(:,i),0.02,'loess');
        r=xcorr(signal(:,i));
        idx=ceil(length(r)/2);
        corr(:,i)=r(idx:end);
        corr_m(i)=mean(r(idx+100:idx+400));
        corr_p(i)=r(idx);
        [peaks,isi]=findpeaks(signal(:,i));
        peaks(peaks<0)=[];
        isi=diff(isi);
        pd=fitdist(peaks,'gamma');
        a(i)=pd.a;
        b(i)=pd.b;
    end
    corr_m=corr_m./max(corr_m);
    
    if nargin==2
        subplot(1,3,1);
    end
    plot(b,corr_m,'.'); hold on;
    xlabel('Gamma b');
    ylabel('Autocorr mean');

    X=[b./max(b);corr_m];
    idx=DBSCAN(X',epsilon,MinPts);

    plot(b(~idx),corr_m(~idx),'x'); hold off
    
    if nargin==2
        subplot(1,3,2);
        imagesc(mask);
        title('Before');
        
        subplot(1,3,3);
        rejects=arrayfun(@(x) find(mask==x),find(~idx),'uniformoutput',false);
%         lol=1:length(mask)^2;
        rejects=cell2mat(rejects)';
        rejected=mask;
        rejected(rejects)=0;
%         lol(rejects)=[];
%         rejects=lol;
        imagesc(rejected);
        title('After');
    end
    
    answer=questdlg('Are you satisfied?','Satisfaction',{'Yes!','NO!'});
    if strcmp(answer,'Yes')
        satisfied=1;
    elseif strcmp(answer,'No')
        thres=inputdlg({'epsilon: ', 'MinPts: '},'DBSCAN Params',1,{num2str(epsilon),num2str(MinPts)});
        epsilon=str2num(thres{1});
        MinPts=str2num(thres{2});
    else
        return
    end
end

tcs.ratio(:,idx==0)=[];
tcs.nCells=size(tcs.ratio,2);
tcs.raw(:,idx==0)=[];
tcs.ring(:,idx==0)=[];
tcs.neuropil(:,idx==0)=[];
tcs.cellIDs(idx==0)=[];
tcs.baseline(idx==0)=[];
tcs.ratio_sd(idx==0)=[];
disp(['Rejected: ' num2str(find(~idx)')]);