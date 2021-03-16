function lp = tz_clustloglk(x,lv,param)
%TZ_CLUSTLOGLK Log likelihood for clusters
%   LP = TZ_CLUSTLOGLK(X,LV) reutrns the log likelihood for clusters. It
%   assumes that each cluster has a normal distribution. The normal
%   distributions are esitimated by MLE. X is a [feature matrix] and lv is
%   a [label vector].
%   
%   LP = TZ_CLUSTLOGLK(X,LV,PARAM) specifies how to estimate the
%   likelihood. PARAM has the following parameters:
%       'model' -
%           'gmx' : gausssin mixture
%               'tz_estnormloglk' - parameters for TZ_ESTNORMLOGLK.
%   
%   See also

%   06-Jun-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('model','gmx','tz_estnormloglk', ...
    struct('covmeth','mlk')));

nsample = size(x,1);
ncluster = max(lv);

lp = 0;
for i=1:ncluster
    clusterSize = sum(lv==i);
    lp = lp+tz_estnormloglk(x(lv==i,:),param.tz_estnormloglk)+ ...
        clusterSize*log(clusterSize);
end

lp = lp-nsample*log(nsample);

