function q=bin_poisson_sim(varargin)
  % Generate simulated Q-matrix modeled as independent Poisson processes
  % Inputs:
  %   m: number of samples
  %   n: number of neurons
  %   'assemblies': logical index of neuronal assemblies with R assemblies
  %       for C neurons
  %   'R': mean spiking probability vector of length n
  %   'prob': vector of individual ensemble activation probability (single value
  %       interpreted as same probability across all ensembles)
  % Outputs:
  %   q: Q-matrix

  m=varargin{1};
  n=varargin{2};
  [assemblies,R,prob]=parse_input(varargin);
  if size(R,1)>1
    R=R';
  end

  q=rand(m,n);
  q=bsxfun(@lt,q,R);

  assemblies=logical(assemblies);

  for i=1:size(assemblies,1)
      events=randi(m,1,round(m*prob(i)));
      q(events,assemblies(i,:))=1;
  end

  function [assemblies,R,prob]=parse_input(inputs)
  assemblies=[];
  R=[];
  prob=[];

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
          otherwise
      end
      idx=idx+1;
  end

  if any([isempty(assemblies) isempty(R) isempty(prob)])
      error('missing inputs');
  end

  if length(prob)==1
    prob=prob.*ones(1,length(R));
  end
