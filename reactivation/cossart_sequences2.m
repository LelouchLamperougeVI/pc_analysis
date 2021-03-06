function sequences=cossart_sequences2(raw,e,assemblies)
% Sequence detection from calcium signals adapted from Malvache et al. 2016
% e of length t with e>0 as assembly number

sampling_rate=19.1; %fps
window=1; %seconds
lag_window=0.5; %seconds
w=round(window*sampling_rate/2);
e_window=round(lag_window*sampling_rate/2);

e=e(e>0);
idx=find(e);

sequences=zeros(length(idx),size(raw,2));

for i=1:length(idx)
  try
    M=raw(idx(i)-w:idx(i)+w,:);
  catch
    try
      M=raw(1:idx(i)+w,:);
    catch
      M=raw(idx(i)-w:end,:);
    end
  end

  reference=median(M(:,assemblies{e(i)}),2);
  cov_map=arrayfun(@(x) xcov(M(:,x),reference,e_window),1:size(M,2),'uniformoutput',false);
  cov_map=cell2mat(cov_map);
  p=arrayfun(@(x) polyfit(1:size(cov_map,1),cov_map(:,x)',2)',1:size(cov_map,2),'UniformOutput',false);
  p=cell2mat(p);
  vertex=-p(2,:)./(2.*p(1,:));
  
  temp=ones(1,length(vertex));
  temp(assemblies{e(i)})=0;
  vertex(logical(temp))=nan;
  
  sequences(i,:)=vertex;
end
