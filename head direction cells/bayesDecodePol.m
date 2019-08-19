function [decoded,Pts]=bayesDecodePol(t,signal,theta_bins)
% Modified bayesian inference decoding for calcium signals with polar behavior

edges=0:2*pi/theta_bins:2*pi;
t_bins=0:2*pi/theta_bins:2*pi-2*pi/theta_bins;
Pt=arrayfun(@(x) length(find(t>=edges(x) & t<edges(x+1))),1:theta_bins);
Pt=Pt./sum(Pt);

if any(~Pt)
    error('Too many bins');
end

% Pt_fun=@(t) Pt(ceil(t./2*pi.*theta_bins));

s=arrayfun(@(x) signal(signal(:,x)>0,x),1:size(signal,2),'uniformoutput',false);
Ps=cellfun(@(x) gamfit(x),s,'uniformoutput',false);
Ps=reshape(cell2mat(Ps),2,size(signal,2));
Ps0=sum(~signal)./size(signal,1);

% Ps_fun=@(s,n) gampdf(s,Ps(1,n),Ps(2,n));


idx=arrayfun(@(x) (t>=edges(x) & t<edges(x+1)),1:theta_bins,'uniformoutput',false);
Pst=cell(size(signal,2),length(idx));
Pst0=zeros(size(signal,2),length(idx));
for i=1:size(signal,2)
    for j=1:length(idx)
        Pst{i,j}=signal(idx{j}' & signal(:,i)>0,i);
        Pst0(i,j)=sum(idx{j}' & signal(:,i)==0)./sum(idx{j});
    end
end
Pst=cellfun(@(x) gamfit(x),Pst,'uniformoutput',false);
Pst=reshape(cell2mat(Pst),size(Pst,1),size(Pst,2),2);

% Pst_fun=@(s,n,t) gampdf(s,Pst(n,t,1),Pst(n,t,2));

Pts=zeros(size(signal,1),size(signal,2),theta_bins);
for i=1:size(signal,1)
    progressbar((i-1)/size(signal,1));
    for j=1:size(signal,2)
        k=1:theta_bins;
        Pts(i,j,k)=Pst_fun(signal(i,j),j,k).*Pt(k)./Ps_fun(signal(i,j),j);
    end
end
progressbar(1);

% Pts(isnan(Pts) | Pts==inf | Pts==0)=1;
Pts=reshape(prod(Pts,2),size(Pts,1),size(Pts,3));

[~,decoded]=max(Pts,[],2);


    function p=Ps_fun(s,n)
        if s>0
            p=gampdf(s,Ps(1,n),Ps(2,n)).*(1-Ps0(n));
        else
            p=Ps0(n);
        end
    end
    function p=Pst_fun(s,n,t)
        if s>0
            p=gampdf(s,Pst(n,t,1),Pst(n,t,2)).*(1-Pst0(n,t));
        else
            p=Pst0(n,t);
        end
    end


end
