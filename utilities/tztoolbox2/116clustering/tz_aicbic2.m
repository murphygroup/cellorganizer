function [aic bic] = tz_aicbic2(feat, cluster, dist, covtype, ispooled)
% TZ_AICBIC2 calculate the AIC/BIC of a clustering result.
%
% [aic bic] = tz_aicbic2(feat, cluster, dist) calculate the AIC/BIC of a
% clustering result. FEAT is a feature matrix, each row of which is
% an observation. CLUSTER could be a cell array, each elementer of which
% gives the members of a cluster by an array of row indices in
% FEAT. CLUSTER could be a [label vector] too. Its elemements correpond to 
% the cluster labels of FEAT. AIC and BIC will be the information criteria 
% calculated from the likelihood and penalty. DIST is the distance function
% for clustering. 
% DIST can be one of the following:
%      'eu' : Euclidean distance (spherical covariance matrix)
%      'ma' : Mahalanobis distance (full covariance matrix)
%
% [aic bic] = tz_aicbic2(feat, cluster) will use the Euclidean
% distance as default, which is the same as [aic bic] =
% tz_aicbic(feat, cluster, 'eu').
%
% [aic bic] = tz_aicbic2(feat, cluster, dist, covtype) specifies the covariance
% type of calculating likelihood by COVTYPE:
%     1 : full covariance matrix
%     2 : diaganol covariance matrix
%     3 : spherical covariance matrix
%     0 : the same as tz_aicbic2(feat, cluster, dist)
%
% [aic bic] = tz_aicbic2(feat, cluster, dist, covtype,1) will use pooled 
% covariance matrix and tz_aicbic2(feat, cluster, dist, covtype,0) is for 
% unpooled. 
%
% See also ml_estpdf, ml_pdf, ml_estcov.

%   Initial written by J. Hua
%   Copyright (c) Center for Bioimage Informatics, CMU
%   $ Date: 2006/06/16 $

% distance function
if ~exist('dist','var')
    dist = 'eu';
end

%initialize parameters for covariace estimation
initm.method = 'fix';

if ~exist('covtype','var')
    covtype = 0;
end

if ~exist('ispooled','var')
    ispooled = 0;
end

if ispooled~=0 & ispooled~=1
    error('ispooled should be either 0 or 1');
end

initm.ispooled = ispooled;

switch(covtype)
    case 1
        initm.tz_estcov.method = 'mle';
    case 2
        initm.tz_estcov.method = 'ind';
    case 3
        initm.tz_estcov.method = 'idv';     
    case 0
        if strcmp(dist(1:2),'eu')
            initm.tz_estcov.method = 'idv';
        elseif strcmp(dist(1:2),'ma')
            initm.tz_estcov.method = 'mle';
        else
        error(['Wrong type of distance function : ' dist]);
        end
    otherwise
        error(['Unrecognized covariace matrix type: ' num2str(covtype)]);
end

[n, nfeat] = size(feat);

if iscell(cluster)
    clusterid = [];
    for i = 1:length(cluster)
        clusterid(cluster{i})=i;
    end
    ncluster = length(cluster);
    if length(clusterid)~=size(feat,1)
        error(['Number of observations not match in feature matrix and' ...
            ' cluster membership cell array.']);
    end
else
    clusterid = cluster';
    ncluster = max(cluster);
end

initm.labels = clusterid';

% if ~ispooled
    % train the stat model
f = ml_gmmfit(feat,initm);
    % calculate the likelihood
    
% else
%     removeidx = [];
%     trainfeat = feat;
%     mu = [];
%     ps = [];
%     for i=1:ncluster
%         idx = find(clusterid==i);
%         if length(idx)==1
%             removeidx = [removeidx;idx];
%         end
%         if isempty(idx)
%             mu(i,:) = zeros(1,size(feat,2));
%         else
%             mu(i,:) = mean(feat(idx,:),1);
%         end
%         
%         trainfeat(idx,:) = ml_addrow(feat(idx,:),-mu(i,:));
%         ps(i) = length(idx)/length(clusterid);
%     end
%     trainfeat(removeidx,:) = [];
% %     f = ml_estpdf(trainfeat, ...
% %         struct('name','mvn'),struct('tz_estcov',initm.tz_estcov));
% %     fsigma = ml_estcov(trainfeat,initm.tz_estcov);
%     [fsigma,mu,ps] = ml_estpoolcov(feat,clusterid);
%     f = struct('name','mix');
%     for i=1:ncluster
%         idx = find(clusterid==i);
%         f.pdfs{i} = struct('name','mvn','mu',mu(i,:),'sigma',fsigma);
%         f.ps = ps;
%     end
    
% end

l = ml_loglk(feat,f);

% calculate the penalty
% penalty for membership


if ~ispooled
    ncov = ncluster;
    p = ncluster-1; 
else
    ncov = 1;
    p = ncluster-1;
end

switch initm.tz_estcov.method
    case 'ind'
    % penalty for mean/variance of each cluster
        p = p+ncluster*nfeat+ncov*nfeat;
    case 'mle'
        p = p+ncluster*nfeat+ncov*nfeat*(nfeat+1)/2;
    case 'idv'
        p = p+ncluster*nfeat+ncov;
end
% calculate the IC
aic = 2*p-2*l;
bic = p*log(n)-2*l;
