function [centers,options,post,errlog] = ...
    tz_kmeans(centers,data,options,removesmall)
%TZ_KMEANS K-means clustering. 
%   CENTERS = TZ_KMEANS(CENTERS,DATA,OPTIONS) does the same thing as
%   KMEANS in NETLAB toolbox does if centers are matrix of seeding
%   coordinates. CENTERS can also be a string of number, like '10'.
%   This specifies the number of clusters and the seeds will be 
%   initialized automatically. To specify seeds initialization method,
%   CENTERS should be a cell array of two strings. The first string is
%   for cluster number as described. And the second string is the option
%   for seeding metod:
%       'rnd' - random
%       'opt' - optimal seeding
%   OPTIONS contains parameters for k-means algorithm. See KMEANS for 
%   details about OPTIONS.   
%
%   CENTERS = TZ_KMEANS(CENTERS,DATA,OPTIONS,REMOVESMALL) will avoid
%   clusters with size no greater than REMOVESMALL. If REMOVESMALL
%   is -1, the size of each cluster will be greater than the number
%   of variables, which is the number of columns of DATA.
%
%   TZ_KMEANS(CENTERS,DATA) does the same thing as
%   KMEANS(CENTERS,DATA,ZEROS(1,14)) does.
%  
%   [CENTERS,OPTIONS,POST,ERRLOG] = TZ_KMEANS(...) also returns other
%   outputs like KEMANS does.

%   18-MAY-2004 Modified T. Zhao
%       - add automatic thresholding
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('At least 2 arguments are required')
end

if ~exist('options')
    options = zeros(1,14);
end

if ~exist('removesmall','var')
    removesmall=0;
end

if iscell(centers)
    seeding = centers{2};
    centers = centers{1};
end

if ~exist('seeding','var')
    seeding = 'opt';
end

if isstr(centers)
    k=str2num(centers);
    n=size(data,1);
    idx = randperm(n);
    switch seeding
    case 'opt'
        idx2(1) = idx(1);
        for m = 2 : k
            dist = [];
            for o = 1 : length(idx2)
                dist(:, o) = sqrt(sum((data - repmat(data(idx2(o), :), n, 1)).^2, 2));
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
    case 'rnd'
        idx2=idx(1:k);
    end
    centers = data(idx2, :);
end

if removesmall==-1
    removesmall=size(data,2);
end

[ndata, data_dim] = size(data);
[ncentres, dim] = size(centers);
% Matrix to make unit vectors easy to construct
id = eye(ncentres);

% Check if centres and posteriors need to be initialised from data
if (options(5) == 1)
    % Assign each point to a randomly-chosen cluster
    index = fix(rand(ndata,1)*ncentres)+1;
    post=id(index,:);
    
    num_points = sum(post, 1)
    % Find the centres based on the random assignments
    for j = 1:ncentres
        if (num_points(j) > 0)
            centres(j,:) = sum(data(find(post(:,j)),:), 1)/num_points(j);
        end
    end
end
%     size(data)
%     size(centers)
[centers,options,post,errlog] = netlab_kmeans(centers,data,options);

if size(data,1)<=removesmall
    warning('the data size is too small');
    return;
end

clustersize=sum(post,1);

while any(clustersize<=removesmall)
    [minsize,idx]=min(clustersize);
    centers(idx,:)=[];   

    [centers,options,post,errlog] = netlab_kmeans(centers,data,options);
    clustersize=sum(post,1);
end

