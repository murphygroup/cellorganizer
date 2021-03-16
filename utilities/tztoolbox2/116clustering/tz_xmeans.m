function [y,centers,bic] = tz_xmeans(x,param)
%TZ_XMEANS X-means clustering.
%   Y = TZ_XMEANS(X) returns a [label vector] of the [feature matrix] X.
%   The labels are memberships of the clusters from x-means algorithm.
%   
%   Y = TZ_XMEANS(X,PARAM) allows setting parameters of clustering. PARAM
%   is a field of the following fields according to xmeans.c writtern by
%   Rudolph van der Merwe. (Details are not included here. You can visit
%   http://www.cs.cmu.edu/~dpelleg/kmeans/readme.txt for details):
%       'k' - Default 1
%       'max_leaf_size' -  Default 40
%       'min_box_width' - Default 0.05
%       'cutoff_factor' - Default 0.5
%       'max_iter' - Default 200
%       'num_splits' - Default 6
%       'max_ctrs' - Default 40
%   
%   [Y,CENTERS,BIC] = TZ_XMEANS(...) also returns cetners and bic value.
%
%   Reference: http://www.cs.cmu.edu/~dpelleg/kmeans/
%
%   See also TZ_CLUSTERING

%   03-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('k',1,'max_leaf_size',10, ...
    'min_box_width',0.05,'cutoff_factor',0.5,'max_iter',200, ...
    'num_splits',6,'max_ctrs',40));

[centers, n, bic] = xmeans(x',param.k,param.max_leaf_size, ...
    param.min_box_width,param.cutoff_factor,param.max_iter, ...
    param.num_splits,param.max_ctrs);

centers = centers';

y = tz_testncclassif(x,struct('centers',centers));
