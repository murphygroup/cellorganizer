function out = tz_kmeans3(x,param)
%TZ_KMEANS3 K-means clustering.
%   OUT = TZ_KMEANS3(X) clusters the [feature matrix] X into 3 clusters 
%   using k-means algorithm. It returns a structure that that contains
%   clustering results:
%       'centers' - centers of the clusters. Each row is a center
%       'options' - the optional parameters. See KMEANS for details.
%       'label' - a vectors of labels. The ith element is the label of the
%           ith row in X.
%       'errlog' - the log of error values. See KMEANS for details.
%       'distance' - distance function. 'euc' or 'mah'.
%   
%   OUT = TZ_KMEANS3(X,PARAM) allows more flexible clustering by specifying
%   PARAM, which is a structure and has following fields:
%       'k' - number of clusters
%       'seeds' - initialization seeds or the way of initianlizing seeds. 
%           If it is a matrix, each row is a seed for one cluster. If it is
%           a string, there are three possibilities. 'sample' selects seeds 
%           randomly. 'uniform' selects seeds uniformly. 'first' takes the
%           first PARAM.k samples as seeds. THe default value is 'sample'.
%       'distance' - distance function. There are two choices: 'euc' for
%           Euclidean distance and 'mah' Mahalanobis distance.
%       'options' - the optional parameters. See KMEANS for details. 
%       'minclustsize' - minimal cluster size.
%   
%   Notice: This function calls KMEANS from the netlab toolbox. Since the
%   statistics toolbox of Matlab 7 also provides KMEANS function, there is
%   a potential conflict of these two funcion and it might cause errors.
%   
%   See also TZ_CLUSTERING

%   24-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('k',3,'seeds','sample','distance','euc', ...
    'options',zeros(1,14),'minclustsize',0));

if param.k==1
    out = struct('centers',mean(x,1), ...
        'label',ones(size(x,1),1), ...
         'dist',param.distance);
    return;
end

nsample=size(x,1);

%Preprocess data for mahalanobis distance
if strcmp(param.distance,'mah')
    covx = cov(x);
    L = sqrtm(inv(covx));
    x = ml_addrow(x,-mean(x,1))*L';
end

%Initialize random seeds
if isstr(param.seeds)    
    switch param.seeds
        case 'uniform'
            idx = randperm(nsample);
            %This part is modified from ML_KMEANS_EUC
            idx2(1) = idx(1);
            for m = 2 :param.k
                dist = [];
              
                for o = 1 : length(idx2)
                    dist(:, o) = sqrt(...
                        sum((x - repmat(x(idx2(o), :),nsample,1)).^2,2));
                end

                [v, I] = sort(sum(dist, 2));
                done = 0;
                c = size(dist, 1);
                while (~done)
                    if (find(idx2 == I(c)))
                        c = c - 1;
                    else
                        idx2(m) = I(c);
                        done = 1;
                    end
                end
            end

        case 'sample'
            idx = randperm(nsample);
            idx2=idx(1:param.k);
        case 'first'
            idx2 = 1:param.k;
    end
    centers = x(idx2, :);
else
    centers = param.seeds;
end

if param.k ~= size(centers,1)
    error('The number of seeds must equal to k');
end

[ndata, data_dim] = size(x);
[ncentres, dim] = size(centers);

[centers,options,post,errlog] = netlab_kmeans(centers,x,param.options);

clustsize = sum(post,1);
while any(clustsize<param.minclustsize)
    [minsize,idx]=min(clustsize);
    centers(idx,:)=[];   
    [centers,options,post,errlog] = netlab_kmeans(centers,x,param.options);
    clustsize=sum(post,1);
end

out = struct('centers',centers,'options',options, ...
    'label',ml_post2label(post),'errlog',errlog, ...
    'dist',param.distance);

