function [features,classlabels,antibodyids] = meddata( subfeatures, sublabels)

data = subfeatures;

antibodyids = unique(sublabels.antibodyids);
features = zeros(length(antibodyids), size(data,2));

idx_ind = zeros(size(antibodyids));
for i=1:length(antibodyids)
    idx = find(antibodyids(i)== sublabels.antibodyids);
    med = median(data(idx,:));
    [a ind] = min(sqrt(sum((data(idx,:)-repmat(med,[length(idx) 1])).^2,2)),[],1);

    idx_ind(i) = idx(ind);
end

features = data(idx_ind,:);
classlabels = sublabels.classlabels(idx_ind);

