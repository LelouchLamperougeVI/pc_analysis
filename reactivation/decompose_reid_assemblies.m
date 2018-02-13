function assemblies=decompose_reid_assemblies(M,assemblies,plotFlag)

  for i=1:length(assemblies)
    temp=M(assemblies{i},:);
    temp=sum(temp);
    idx=kmeans(temp',2);
    idx=logical(idx-1);
    if mean(temp(idx))<mean(temp(~idx))
      idx=~idx;
    end
    assemblies{i}=find(idx);
  end

  if plotFlag
      fprintf('Detected assemblies: \n');
      arrayfun(@(x) fprintf(['\t Assembly ' num2str(x) ':\t' mat2str(assemblies{x}') '\n']),1:length(assemblies));
  end
