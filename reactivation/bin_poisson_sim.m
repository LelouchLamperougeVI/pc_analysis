function [q,e]=bin_poisson_sim(varargin)
% Generate simulated Q-matrix modeled as independent Poisson processes
% Inputs:
%   m: number of samples
%   n: number of neurons
%   'assemblies': logical index of neuronal assemblies with R assemblies
%       for C neurons
%   'R': mean spiking probability vector of length n
%   'prob': vector of individual ensemble activation probability (single value
%       interpreted as same probability across all ensembles)
%   'miss': mean number of missing members during ensemble activation
%   'jitter': s.d. of normal probability distribution of temporal jitters
%   'jitter_type': type of jittering noise ('normal' - default or
%                  'uniform')
%   'corr': degree of correlation between ensemble members
%   'seed': deconv used for seeding calcium spikes
% Outputs:
%   q: Q-matrix
%   e: reactivation times

m=varargin{1};
n=varargin{2};
[assemblies,R,prob,miss,jitter,jitter_type,C,seed]=parse_input(varargin);
if size(R,1)>1
    R=R';
end

if isempty(C)
    dff_train=1;
else
    if isempty(seed)
        dff_train=randn(m,n)+10;
    else
        sm=min(seed(seed~=0));
        seed=log(seed);
        seed(isinf(seed))=nan;
        seed=(seed-mean(seed,'omitnan'))./std(seed,'omitnan');
        seed=seed(~isnan(seed));
        seed=exp(seed);
        seed=seed-min(seed)+sm;
        dff_train=seed(randi(length(seed),1,m*n));
        dff_train=reshape(dff_train,m,n);
    end
%     C=(assemblies'.*C)*(assemblies.*C);
%     temp=assemblies'*assemblies;
%     C=temp./max(temp(:)).*C;
%     for i=1:n; C(i,i)=1; end
    temp=assemblies'*assemblies;
    temp=temp + diag((max(temp(:)) - max(temp(:)).*C)./C .*ones(1,length(temp)));
    C=diag(1./sqrt(diag(temp))) * temp * diag(1./sqrt(diag(temp)));
    L=chol(C);
    dff_train=dff_train*L;
end

q=rand(m,n);
q=bsxfun(@lt,q,R);
e=zeros(m,size(assemblies,1));
if ~isempty(assemblies) && ~isempty(prob)
    assemblies=logical(assemblies);
    
    for i=1:size(assemblies,1)
        events=nr_randi(m,1,round(m*prob(i)));
        for j=1:length(events)
            nomembers=poissrnd(miss(i));
            idx=find(assemblies(i,:));
            if nomembers>length(idx)
                nomembers=length(idx);
            end
            idx(nr_randi(length(idx),nomembers,1))=[];
            
            try
                if strcmpi(jitter_type,'normal')
                    shift=randn(1,length(idx));
                    shift=floor(abs(shift).*jitter).*(abs(shift)./shift);
                elseif strcmpi(jitter_type,'uniform')
                    shift=randi(2*jitter+1,1,length(idx))-jitter-1;
                else
                    error('undefined jitter_type');
                end
                shift=sub2ind([m n],shift+events(j),idx);
            catch
                shift=sub2ind([m n],ones(1,length(idx)).*events(j),idx);
            end
            
            q(shift)=1;
        end
        e(events,i)=1;
    end
end
q=q.*dff_train;


function [assemblies,R,prob,miss,jitter,jitter_type,C,seed]=parse_input(inputs)
assemblies=[];
R=[];
prob=[];
miss=0;
jitter=0;
jitter_type='normal';
C=[];
seed=[];

idx=3;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'assemblies'
            idx=idx+1;
            assemblies=inputs{idx};
        case 'r'
            idx=idx+1;
            R=inputs{idx};
        case 'prob'
            idx=idx+1;
            prob=inputs{idx};
        case 'miss'
            idx=idx+1;
            miss=inputs{idx};
        case 'jitter'
            idx=idx+1;
            jitter=inputs{idx};
        case 'jitter_type'
            idx=idx+1;
            jitter_type=inputs{idx};
        case 'corr'
            idx=idx+1;
            C=inputs{idx};
        case 'seed'
            idx=idx+1;
            seed=inputs{idx};
        otherwise
    end
    idx=idx+1;
end

if isempty(R)
    error('missing spiking probability input');
end

if length(prob)==1
    prob=prob.*ones(1,size(assemblies,1));
end

if length(R)==1
    R=R.*ones(1,inputs{2});
end

if length(miss)==1
    miss=miss.*ones(1,size(assemblies,1));
end
