function labels = tz_kmeans_mah(centers,feat,param)
%TZ_KMEANS_MAH K-means clustering using Mahalanobis distance.
%   LABELS = TZ_KMEANS_MAH(CENTERS,FEAT) returns the clustering [label 
%   vector] of the [feature matrix] FEAT. CENTERS is a [feature matrix] for
%   initial centers. The number of clusters is the number of rows of
%   CENTERS.
%   
%   LABELS = TZ_KMEANS_MAH(CENTERS,FEAT,PARAM) allows specifying parameters
%   for clustering:
%       'maxiter' - maximum iteration (default 100)
%       'e' - error for stop (default 1e-5)
%
%   See also

%   23-Mar-2006 Modified from ML_KMEANS_MAH T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('maxiter',100,'e',1e-5));

%No. of instances
n = size(feat, 1);
k = size(centers,1);

clust = zeros(n, k);
id = eye(k);
% clust(idx2,:) = id;
done = 0;
iter = 0;
while (~done)
    iter = iter+1;
    dist = [feat;centers];
    cen2 = centers;

    pd = squareform(pdist(dist, 'mahal'));
    pd = pd(n + 1 : n + k, 1 : n);

    [Y, I] = min(pd, [], 1);
    clust = id(I, :);
    no_instances = sum(clust, 1);

    for m = 1 : k
        if (no_instances(m) > 0)
            centers(m, :) = mean(feat(find(clust(:,m)),:), 1);
        end
    end
    
    if iter>param.maxiter
        warning('maximum iteration exceeded');
        done = 1;
    end

    if max(abs(centers(:)-cen2(:)))<=param.e
        done = 1;
    end
end

for m = 1 : k
    tmp = clust(:, m);
    labels(find(tmp),:) = m;
end
