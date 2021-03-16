function [ output_args ] = proj_error_eig2( model )

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

niter = 50;

mse = cell(1,niter);


figure('color', 'w'), 
for n = 1:niter

    inds = randperm(size(dists,1));
    
    pos_prev = eig(dists(inds(1:2), inds(1:2)), 'vector');

    for i = 3:length(dists)
        pos = eig(dists(inds(1:i), inds(1:i)), 'vector');
        pos = pos / pos(1);


        [~, ind] = min(abs(pos));

        pos_prev = [pos_prev(1:ind-1); 0; pos_prev(ind:end)];
        pos_prev = pos_prev / pos_prev(1);

        mse{n}(i-2) = norm(pos - pos_prev,2);

        pos_prev = pos;
    end

    plot(mse{n}, 'k')
    hold on
end

mu = mean(vertcat(mse{:}));
sigma = std(vertcat(mse{:}));
plot(mu, 'r')
plot(mu+sigma, '--r')
plot(mu-sigma, '--r')

xlabel('Number of shapes in complete shape space')
ylabel('error from previous iteration')


end

