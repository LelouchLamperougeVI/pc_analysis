function corr(obj)
% build correlation matrix

deconv = obj.twop.deconv;
% deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
% deconv(:, all(isnan(deconv))) = 0;
deconv=fast_smooth(deconv,obj.ops.sig*obj.twop.fs);
% deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
deconv = deconv(~isnan(deconv));
deconv = reshape(deconv, [], size(obj.twop.deconv, 2));
deconv=(deconv-mean(deconv))./std(deconv); %zscore

obj.ensembles.R=corr(deconv);

% sig=obj.ops.sig*obj.twop.fs;
% shuffles = obj.ops.shuffles;

% deconv = obj.twop.deconv;
% deconv(isnan(sum(deconv,2)),:)=[];
% null_R = zeros(size(deconv,2),size(deconv,2),obj.ops.shuffles);

% dq=parallel.pool.DataQueue;
% afterEach(dq,@updateBar);
% h=waitbar(0,'permutation testing...');
% count=1;
% parfor i=1:obj.ops.shuffles
%     temp=burst_shuffler(deconv);
%     temp=(temp-mean(temp,'omitnan'))./std(temp,'omitnan'); %zscore
%     temp=fast_smooth(temp,sig);
%     null_R(:,:,i) = corr(temp);
%     send(dq,i);
% end
% obj.ensembles.null_R = null_R;
% obj.ensembles.null_R = mean(obj.ensembles.null_R,3);
% close(h);

%     function updateBar(~)
%         waitbar(count/shuffles,h);
%         count=count+1;
%     end
end
%% Old version 
% compute the corr matrix inside and outside SCE events separately
% consider corr matrix obtained outside SCE as null
% doesn't work so good in noisy data

% idx=[];
% for i=1:length(obj.ensembles.SCE.on)
%     idx=[idx find(obj.ensembles.SCE.on(i)==obj.ts) : find(obj.ensembles.SCE.on(i)+obj.ensembles.SCE.dur(i)==obj.ts)];
% end
% 
% deconv=obj.deconv;
% deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
% deconv=fast_smooth(deconv,obj.ops.sig*obj.fs);
% 
% obj.ensembles.R=corr(deconv(idx,:));
% 
% idx=setxor(idx,1:size(deconv,1));
% idx=setxor(idx, find(isnan(sum(deconv,2))));
% 
% obj.ensembles.null_R=corr(deconv(idx,:));