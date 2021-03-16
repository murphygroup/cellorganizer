function [feats,names] = spatialFeatures( img1, img2, img3, DSF)


p = spatialStatistics( img1, img2, DSF); 
nbins = length(p);


feats = p;
% pn = spatialStatistics( img1, img3, DSF); 

% feats = p./pn;
% feats = feats / sum(feats) * (nbins-1);


names = [repmat({'histfeats'},[1 nbins]) ];
b = 0:1/(nbins-1):1;
for i=1:nbins+1
    idx = i:nbins+1:length(feats);
    for j=1:length(idx)
        names{idx(j)} = [names{idx(j)} '_' num2str(b(i))];
    end
end

return