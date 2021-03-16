function [ output_args ] = proj_error_procrustes( model )

dists = model.cellShapeModel.distances_incomplete;
dists(eye(size(dists))>0) = nan;

keepinds = find(~all(isnan(dists),2));

dists = dists(keepinds,keepinds);
dists(eye(size(dists))>0) = 0;

 nancounts = sum(isnan(dists),2);

while any(nancounts)
    [~, ind] = max(nancounts);
    
    dists(ind,:) = [];
    dists(:, ind) = [];
    keepinds(ind) = [];
    
    nancounts = sum(isnan(dists),2);
end

pos = cmdscale(dists);

pos_prev = cmdscale(dists(1:2, 1:2), 1);

for i = 3:length(dists)
    pos = cmdscale(dists(1:i, 1:i), i-1);
    mse(i-2) = procrustes_mse(pos_prev, pos(1:end-1,:));
    
    pos_prev = pos;
end


end

