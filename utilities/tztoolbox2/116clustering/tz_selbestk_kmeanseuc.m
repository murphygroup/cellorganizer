function results=tz_selbestk_kmeanseuc(combobj,ntrial,randseed)
%TZ_SELBESTK_KMEANSEUC Find optimal clusters from k-means.
%   RESULTS = TZ_SELBESTK_KMEANSEUC(X,NTRIAL,RANDSEED) resturns a
%   structure containing the information about the optimal clusters.
%   RANDSEED is the random state for rand function to make results
%   repeatable.
%   RESULTS has the following fields (see KMEANS for more details):
%       k - number of clusters
%       centers - centeral coordinates of each cluster
%       post - cluster membership of each sample, 1-N coding
%       errlog - error
%       aic - aic value
%       clusterinfo - cluster information from aic calculation

%   ??-???-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%normalize them all
fnorm=mb_featurenorm(combobj,zeros(size(combobj)));
clear combobj

% Seed the random number generator
%randseed = [1:35]'+2; %zeros(1,35);
rand('state',randseed);

for tr = 1:ntrial
    if( tr == 1)
        first_k = 1;
    else
        first_k = 1;
    end
    for k = first_k:40
        %% choose a random subset of samples for seeds
        randidx = randperm(size(fnorm,1));
        seeds = fnorm(randidx(1:k),:);
        
        % Do k-means clustering
        options = zeros(1,14);
        %options(5) = 1; % let seeds be picked randomly by rm_mahalkmeans
        %seeds = ones(k,size(fnorm,2));
        [centers,options,post,errlog]=netlab_kmeans(seeds,fnorm,options);
        if( size(post,2) == size(centers,1))
            actual_k = size(post,2);
        else
            error('size(post,2) is not the smae as size(centers,1))');
        end
        
        % Judge goodness of this k by Akaike Information Criterion
        [aic, clusterinfo] = rm_akaike_euclid( fnorm, post);
        if(k==first_k & tr==1)
            minaic=aic;
        end
        
        %AICs = [AICs aic];
        %Ks = [Ks actual_k];
        %figure(1)
        %plot(Ks,AICs);
        %figure(1)
        if aic<minaic
            minaic=aic;
            results = struct('k',actual_k,'centers',centers,'post',post,'errlog',errlog, ...
                'aic',aic,'clusterinfo',clusterinfo);
        end
        
       % save(['clustering_results_euclid/kmeanseuc_aiceuc_f11_tr' ...
       %         num2str(tr,'%.2d') '_k' num2str(actual_k,'%.2d')],'results');
    end
end
