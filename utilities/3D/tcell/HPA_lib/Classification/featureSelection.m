function idx = featureSelection( traindata, trainlabels, method, uniquename)

U = unique( trainlabels);
feat = [];
for j=1:length(U)
    idx = find(trainlabels==U(j));
    feat{j} = traindata(idx,:);
end

switch method
  case 'sda'
    idx = ml_stepdisc( feat, ['idxsda' uniquename '.txt']);
    delete(['idxsda' uniquename '.txt']);
  case 'none'
    idx = 1:size(traindata, 2);
  otherwise
    error( 'this method of feature selection is not supported yet');
end
