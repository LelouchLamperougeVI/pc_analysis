function sequences=pca_sequences(pc,deconv)
% get replay sequences from PC

trials=pca_trials(pc,1);

for i=1:length(trials)
    ranks=[];
    for j=1:size(trials{i},2)
        M=deconv(trials{i}(1,j):trials{i}(2,j),:);
        blanks=any(M); % make sure cell fired in trial
        [~,M]=max(M);
        M=M.*blanks;
        ranks=[ranks;M];
    end
    tmp=double(~logical(ranks));
    ranks=reshape(ranks,1,[]);
    ranks(ranks==0)=[];
    ranks=tiedrank(ranks);
    count=1;
    for h=1:size(tmp,2)
        for k=1:size(tmp,1)
            if tmp(k,h)==1
                tmp(k,h)=ranks(count);
                if h~=size(tmp,2) && k~=size(tmp,1)
                    count=count+1;
                end
            end
        end
    end
    ranks=tmp;
    sequences{i}=sum(ranks)./sum(logical(ranks));
    [~,sequences{i}]=sort(sequences{i});
end