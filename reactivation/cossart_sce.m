function e=cossart_sce(q,thres)

  if nargin<2
    thres=2;
  end

  number_sims=1000;
  C=zeros(1,number_sims*size(q,1));
  for i=1:number_sims
      temp=q;

      shift=randi(size(temp,1),1,size(temp,2));
      temp=mat_circshift(double(temp),shift);

      C((i-1)*size(temp,1)+1:i*size(temp,1))=sum(temp,2);
  end

  [mu,sigma]=normfit(C);
  sd_series=(sum(q,2)-mu)./sigma;
  e=sd_series>thres;
