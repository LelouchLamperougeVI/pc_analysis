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
  % Outputs:
  %   q: Q-matrix
  %   e: reactivation times

  m=varargin{1};
  n=varargin{2};
  [assemblies,R,prob,miss]=parse_input(varargin);
  if size(R,1)>1
    R=R';
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
          idx(nr_randi(length(idx),nomembers,1))=[];
          q(events(j),idx)=1;
        end
        e(events,i)=1;
    end
  end



function [assemblies,R,prob,miss]=parse_input(inputs)
  assemblies=[];
  R=[];
  prob=[];
  miss=0;

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
