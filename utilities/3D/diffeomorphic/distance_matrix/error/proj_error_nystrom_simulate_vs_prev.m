function [ output_args ] = proj_error_nystrom_simulate_vs_prev( model )
%Compares the l2-norm of the eigenvalues for an nystrom embedded matrix 
%over multiple iterations of adding landmarks


dists_complete = model.cellShapeModel.distances_incomplete;
dists_complete(eye(size(dists_complete))>0) = nan;

keepinds = find(~all(isnan(dists_complete),2));

dists_complete = dists_complete(keepinds,keepinds);
dists_complete(eye(size(dists_complete))>0) = 0;

 nancounts = sum(isnan(dists_complete),2);

while any(nancounts)
    [~, ind] = max(nancounts);
    
    dists_complete(ind,:) = [];
    dists_complete(:, ind) = [];
    keepinds(ind) = [];
    
    nancounts = sum(isnan(dists_complete),2);
end

niter = 50;

for j = 1:niter
    j
    inds = randperm(size(dists_complete,1));
    dists_complete = dists_complete(inds,:);
    dists_complete = dists_complete(:,inds);
    
    keepinds = false(length(dists_complete));
    keepinds(1:2) = true;

    dists_temp = nan(size(dists_complete));
    dists_temp(keepinds,:) = dists_complete(keepinds,:);
    dists_temp(:,keepinds) = dists_complete(:,keepinds);

    dists_temp(logical(eye(size(dists_temp)))) = 0;


    [pos_prev, ~, info] = embed_partial_distance_matrix(dists_temp, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', false));

    dists_prev = squareform(pdist(pos_prev));
    
    eigval_prev = sqrt(info.unfiltered_result_eigenvalues);

    for i = 3:length(dists_complete)
        
        keepinds(i) = true;

        dists_temp(keepinds,:) = dists_complete(keepinds,:);
        dists_temp(:,keepinds) = dists_complete(:,keepinds);

        [pos, ~, info] = embed_partial_distance_matrix(dists_temp, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', true, 'desired_shape_space_dimensionality', inf));

        dists = squareform(pdist(pos));
        
        eigval = sqrt(info.unfiltered_result_eigenvalues');

        eigval = eigval / eigval(1);
        [~, ind] = min(abs(eigval));

        d1 = length(eigval);
        d2 = length(eigval_prev);
        delta = d1-d2;

        if d1-d2 < 0
            [~, ind] = min(abs(eigval));
            eigval_temp = [eigval(1:ind-1); zeros(abs(delta),1); eigval(ind:end)];
            eigval_prev_temp = eigval_prev;
        else
            [~, ind] = min(abs(eigval_prev));
            eigval_prev_temp = [eigval_prev(1:ind-1); zeros(abs(delta),1); eigval_prev(ind:end)];
            eigval_temp = eigval;
        end

        eignorm(j,i-2) = norm(eigval_temp - eigval_prev_temp,2);

        matnorm(j,i-2) = norm(dists - dists_prev,2);
        matnorm_fro(j,i-2) = norm(dists - dists_prev, 'fro');
        matmse(j,i-2) = sqrt(sum((dists(:) - dists_prev(:)).^2));
        
        dists_prev = dists;
        eigval_prev = eigval;
    end
end


mu = mean(eignorm,1);
st = std(eignorm,1);


figure('color', 'w')
plot(eignorm', 'k')
hold on,

plot(mu, 'r')
plot(mu+st, '--r')
plot(mu-st, '--r')

axis([0, size(dists_complete,1), 0, prctile(eignorm(:,1), 95)])

end

