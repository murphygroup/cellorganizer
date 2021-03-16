function out = tz_clustering(x,param)
%TZ_CLUSTERING Clustering.
%   OUT = TZ_CLUSTERING(X) returns a structure which represents the output
%   of k-means clustering of the [feature matrix] X. The seeds of the
%   k-means clustering are selected randomly. X will be zscored before
%   clustering. OUT may have different fields in different methods, but it
%   always has a field 'label', which is a column vector of clustering
%   labels of the data. This will output 3 clusters. To customize the
%   number of clusters, please see comments below.
%   
%   OUT = TZ_CLUSTERING(X,PARAM) lets the user specify parameters for
%   clustering. PARAM is a structure and has following fields:
%       'method' - clustering method. Currently it supports following
%           methods:
%               'kmeans' : kmeans algorithm (default).
%               'kmeansaic' : kmeans by aic selection
%               'kmeansbic' : kmeans by bic selection
%               'xmeans' : xmeans algorithm. See XMEANS.
%       'zscore' - 1 for zscoring before clustering (default). otherwise 
%           no zscore will be done.
%       'transform' - transformation of the data.
%           'none': no transformation (default).
%           'pca': pca transformation.
%               'tz_data2pca' - see TZ_DATA2PCA for more details.
%           'spectral': transformation for spectral clustering
%               'ml_spectral' - see ML_SPECTRAL for more details.
%       The following fileds depend on 'method':
%       'kmeans','kmeansaic','kmeansbic'
%          'kmeans' - a structure of kmeans parameters. See TZ_KMEANS3.
%          'k' - number of clusters.
%       'kmeansaic','kmeansbic'
%           'aicbic' - 'new' or 'old'
%           'ks' - numbers of clusters to search. default 2:10.
%           'tr' - number of trials. default 10.
%           'covtype' - covariance type. See TZ_AICBIC2. The default value
%               is 0.
%           'ispooled' - pooled covariance matrix or not. See TZ_AICBIC2.
%               The default value is 0.
%       'xmeans'
%           'xmeans' - a structure of xmeans parameters. See TZ_XMEANS.
%   
%   Example:
%       Cluster data X into K clusters:
%           tz_clustering(X,struct('kmeans',struct('k',K)));
%
%   See also TZ_KMEANS3 TZ_XMEANS

%   24-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU


if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

%Notice: This is coupled with tz_loadclusterdir. Make sure consistent
%updates are made in the future.

param = ml_initparam(param, ...
    struct('method','kmeans','zscore',1,'transform','none', ...
    'kmeans',struct([]),'ml_spectral',struct([]),'tz_data2pca',struct([])));

if any(isnan(x(:)) | isinf(x(:)))
    error('The input data can not be clustered because they contain Nan or Inf values.');
end

if param.zscore==1
    [x,out.zmean,out.zsdev]=ml_zscore(x);
%     x = zscore(x);
end

switch param.transform
    case 'pca'
        x = tz_data2pca(x,param.tz_data2pca);
    case 'spectral'
        x = ml_spectral(x,param.ml_spectral);
    case 'none'
        %do nothing
    otherwise
        error('Unrecognized transformation method.');
end

switch param.method
    case 'kmeans'
        out = tz_kmeans3(x,param.kmeans);
    case {'kmeansaic','kmeansbic'}
        param = ml_initparam(param, ...
            struct('ks',2:10,'tr',10, ...
            'aicbic','new','covtype',0,'ispooled',0));
        inc = 1;
        for k=param.ks
            k
            for t=1:param.tr
                if isempty(param.kmeans)
                    param.kmeans = struct('k',k,'seeds','sample');
                else
                    param.kmeans.k = k;
%                     param.kmeans.seeds = 'sample';
                    param.kmeans = ml_initparam(param.kmeans,struct('seeds','sample'));
                end
                param2 = param;
                if strcmp(param.kmeans.seeds,'first')
                    param2.kmeans.seeds = x((1:param.kmeans.k)+t,:);
                end
                if isnumeric(param.kmeans.seeds)
                    if size(param.kmeans.seeds,2)==1
                        seedidx = ...
                            randsample(param.kmeans.seeds,param.kmeans.k);
                        param2.kmeans.seeds = x(seedidx,:);
                    end
                end
                tmpout{t,inc} = tz_kmeans3(x,param2.kmeans);
                
                switch param.aicbic
                    case 'new'
                        [aics(t,inc),bics(t,inc)] = tz_aicbic2(x, ...
                             tmpout{t,inc}.label,tmpout{t,inc}.dist, ...
                             param.covtype,param.ispooled);
                    case 'old'
                        [aics(t,inc),bics(t,inc)] = tz_aicbic(...
                             x,tmpout{t,inc}.label,tmpout{t,inc}.dist);
                    case 'most'
                        [aics(t,inc),bics(t,inc)] = tz_aicbic3(x, ...
                             tmpout{t,inc}.label,tmpout{t,inc}.dist, ...
                             param.covtype,param.ispooled);
                end
            end
            inc = inc+1;
        end
        if strcmp(param.method,'kmeansaic')
%             [minaic,minidx] = min(mean(aics,1));
            [minaic,tridx,kidx] = ml_min(aics);
        else
%             [minbic,minidx] = min(mean(bics,1));
            [minbic,tridx,kidx] = ml_min(aics);
        end
        out.aics = aics;
        out.bics = bics;
        out.bestk = param.ks(kidx);
        out.label = tmpout{tridx,kidx}.label;
        out.alllabels = tmpout;
        out.dist = tmpout{1}.dist;
    case 'xmeans' %x-means algorithm
        if ~isfield(param,'tz_xmeans')
            param.tz_xmeans = struct([]);
        end
        
        [y,centers,bic] = tz_xmeans(x,param.tz_xmeans);
        out.bic = bic;
        out.label = y;
        out.centers = centers;
end
    

