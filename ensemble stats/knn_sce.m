function [sce,cost]=knn_sce(deconv,assemblies)

k=1; %hard limit for now...

% figure;
% plot(deconv(:,i));
% hold on
% plot(deconv(:,j));
sce=cell(1,length(assemblies));
for a=1:length(assemblies)
    dec=deconv(:,assemblies{a});
    dec=ca_filt(dec);
    dec=log(dec);
    dec(isinf(dec))=nan;
    dec=(dec-mean(dec,'omitnan'))./std(dec,'omitnan');
    
    dec=avg_smooth(dec,5); %smooth with averaging kernel to increase local correlation

    %     figure;
    sce{a}=false(size(dec,1),size(dec,2),size(dec,2));
    for i=1:size(dec,2)
        for j=i+1:size(dec,2)
            %         sce=false(size(deconv,1),1);
            pre=[];I=[];
%             count=1;
            converge=false;
            overlap=false(size(dec,1),1);
            while(~converge)
                pre=[pre I];
                
                idx=find(~sce{a}(:,i,j));
                
                sce{a}(idx(overlap),i,j)=true;
                B=dec(~sce{a}(:,i,j),[i j]);
                A=B;
                
                term=sum(sum(isnan(A)).*log(sum(isnan(A))./size(B,1)))-sum(all(isnan(A),2))*log(sum(all(isnan(A),2))/size(B,1));
                
                A=A(all(~isnan(A),2),:); %continuous-continuous space
                if isempty(A)
                    warning('Empty differential joint entropy set detected. Cannot converge current sample.');
                    break
                end
                
                d=knnsearch(A,A,'k',k+1);
                d=d(:,end); % k th neighbor index
                d=A(d,:); %each row's k-nn
                dist=abs(d-A);
                
                d=abs(B-permute(A,[3 2 1]));
                n=sum(d<permute(max(dist,[],2),[2 3 1]),1);
                overlap=sum(d<permute(max(dist,[],2),[2 3 1]),3);
                converge=max(overlap(~isnan(B) & ~all(~isnan(B),2)))<=1;
                try
                    overlap=overlap==max(overlap(~isnan(B) & ~all(~isnan(B),2))) & ~isnan(B) & ~all(~isnan(B),2);
                catch
                    break
                end
                %     overlap=xor(overlap(:,1),overlap(:,2));
                overlap=any(overlap,2);
                %     overlap(all(~isnan(B(:,[i j])),2))=false;
                
                %         overlap=find(overlap,randi(sum(overlap)));
                %         overlap=overlap(end);
                
                n=permute(n,[3 2 1]);
                n=n./size(B,1).*size(A,1);
                k_xyz=mean(sum(psi(n+1),2));
                
                I=(psi(k)-k_xyz+psi(size(A,1)))*(size(A,1)/size(B,1))-term/size(B,1);
                
                if length(pre)>1
                    converge=converge || (I-pre(end))<median(diff(pre));
                end
                
%                                 plot(count,I,'o'); count=count+1;
                %                 hold on
                
                %     plot(sce.*(10+count*2)); count=count+1;
            end
            sce{a}(:,j,i)=sce{a}(:,i,j);
            %             hold off
        end
    end
    sce{a}=sce{a}.*~isnan(dec);
    sce{a}=sum(sce{a},3);
    sce{a}=sce{a}>arrayfun(@(x) prctile(sce{a}(sce{a}(:,x)>0,x),15),1:size(sce{a},2));
    sce{a}=~sce{a};
    
    disp(['Ensemble (' num2str(a) '/' num2str(length(assemblies)) ') completed']);
end