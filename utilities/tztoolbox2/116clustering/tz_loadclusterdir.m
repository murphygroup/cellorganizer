function [aics,bics,nclust] = tz_loadclusterdir(resultdir,icupdate)
%TZ_LOADCLUSTERDIR Load clustering files from a directory.
%   AICS = TZ_LOADCLUSTERDIR(RESULTDIR) read mat files under the directory
%   RESULTDIR and get clustering information. The file should be named with
%   'cluster(k)', where k is the number of clusters. It will return aic
%   values of the clusters.
%   
%   AICS = TZ_LOADCLUSTERDIR(RESULTDIR,ICUPDATE) will recaluate AICS base
%   on the parameters specified by ICUPDATE:
%       'method' - 'old' (TZ_AICBIC) or 'new' (TZ_AICBIC2)
%           'old': TZ_AICBIC will be called
%               'distfun' - distance function (See TZ_AICBIC)
%           'new': TZ_AICBIC2 will be called
%               'dist' - distance function default 'eu'
%               'covtype' - covariance type. default 0.
%               'pool' - pooled corvariance or not. default 0.
%       'features' - a [feature matrix]  for calculating ICs
%            If 'features' does not exist, the features will be loaded 
%            automtically from the feature file. 
%       'recalc' -  if it is 0, the function will not do the calculation if 
%            the same updating has been done before. If it is 1, it will
%            calculate it anyway. If it is -1, the correponding clusters
%            will be skipped if the update does not exsit. The default
%            value is 0.
%       'preprocess' -  if it is 0, the function will calculate the aics and
%            bics based on the original input features. If it is 1, the 
%            features will be preprocessed (to the 'tranform' stage) before
%            IC calculation. This is only useful when 'features' is specified. 
%            The default value is 0.
%       NOTICE: Be sure to set preprocess to 1 when 'features' is specified.
%            Or you will get wrong results. The default value of 'preprocess'
%            is not set to 1 because of some historical reasons.
%            
%   [AICS,BICS,NCLUST] = TZ_LOADCLUSTERDIR(...) also return bic values and
%   number of clusters.
%   
%   See also

%   18-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

if ~exist('icupdate','var')
    icupdate = struct([]);
end

%icupdate = ml_initparam(icupdate,struct('recalc',0,'preprocess',0));

files = tz_sortfile(tz_ls([resultdir filesep 'cluster*.mat']));

if isempty(files)
    error(['The directory ' resultdir ' has no mat file']);
else
    if ~isempty(icupdate)    
        if ~isfield(icupdate,'features')
            icupdate.preprocess = 1;
        end
        
        if icupdate.preprocess==1
            tmp = load([resultdir filesep files{1}]);
            tokens = tz_strtok(tmp.script,'''');
            featureFile = tokens{2};
            paramFile = tokens{4};
            if ~isfield(icupdate,'features')
                tmp=load(featureFile);
                icupdate.features = tmp.features;
            end
            tmp=load(paramFile);
            tmp.clustering = ml_initparam(tmp.clustering, ...
                struct('zscore',1,'transform','none', ...
                'ml_spectral',struct([]),'tz_data2pca',struct([])));
        end
        if icupdate.preprocess==1           
            if tmp.clustering.zscore==1
                icupdate.features = ml_zscore(icupdate.features);
            end
            
            switch tmp.clustering.transform
                case 'pca'
                    icupdate.features = ...
                        tz_data2pca(icupdate.features, ...
                        tmp.clustering.tz_data2pca);
                case 'spectral'
                    icupdate.features = ml_spectral( ...
                        icupdate.features,tmp.clustering.ml_spectral);
            case 'none'
                %do nothing
            end
        end
    end
end
aics = [];
bics = [];
nclust = [];
for i=1:length(files)
    result = load([resultdir filesep files{i}]);
    i
    if ~isempty(icupdate)
        icupdate = ml_initparam(icupdate,struct('recalc',0,'preprocess',0));
        iscalc = 1; %recalculate aic,bic or not
        %initialize parameters
        switch icupdate.method
            case 'old'
                icupdate = ml_initparam(icupdate,struct('distfun', ...
                                                        'euc'));
            case 'new'
                icupdate = ml_initparam(icupdate, ...
                    struct('dist','eu','covtype',0,'pool',0));
            case 'most'
                icupdate = ml_initparam(icupdate, ...
                    struct('dist','eu','covtype',0,'pool',0));
            otherwise
                error(['Unrecognized method: ' icupdate.method]);
        end
        
        param = rmfield(icupdate,{'features','preprocess','recalc'});
        
        matchedIcIndex = 0;
        %check if the specified ics have been calculated
        if isfield(result,'updateics')
            for k=1:length(result.updateics)
                if ml_structcmp(result.updateics{k}.param,param)==1
                    if icupdate.recalc<=0
                        iscalc = 0;
                    end                   
                    matchedIcIndex = k;
                    break
                end
            end
        end
        
        if matchedIcIndex==0
            if icupdate.recalc==-1
                continue;
            end
        end

        if iscalc==1 %recalculate
            clear allaics allbics
            switch icupdate.method
                case 'old'
                    for j=1:length(result.cluster.alllabels)
                        [allaics(j),allbics(j)] = ...
                            tz_aicbic(icupdate.features, ...
                            result.cluster.alllabels{j}.label, ...
                            icupdate.distfun);
                        %max(result.cluster.alllabels{j}.label(:))
                    end
                case 'new'
                    for j=1:length(result.cluster.alllabels)
                        [allaics(j),allbics(j)] = ...
                            tz_aicbic2(icupdate.features, ...
                            result.cluster.alllabels{j}.label, ...
                            icupdate.dist,icupdate.covtype,icupdate.pool);
                    end
                case 'most'

                    for j=1:length(result.cluster.alllabels)
                        [allaics(j),allbics(j)] = ...
                            tz_aicbic3(icupdate.features, ...
                            result.cluster.alllabels{j}.label, ...
                            icupdate.dist,icupdate.covtype,icupdate.pool);
                    end

                otherwise
                    error(['Unrecognized method: ' icupdate.method]);
            end
            if isfield(result,'updateics')
                if matchedIcIndex==0
                    result.updateics{end+1} = struct(...
                        'param',param,'aics',allaics,'bics',allbics);
                else
                    result.updateics{matchedIcIndex} = struct(...
                        'param',param,'aics',allaics,'bics',allbics);
                end
            else
                result.updateics = ...
                    {struct('param',param,'aics',allaics,'bics',allbics)};
            end
            
            ml_savestruct([resultdir filesep files{i}],result);
        else
            allaics = result.updateics{matchedIcIndex}.aics;
            allbics = result.updateics{matchedIcIndex}.bics;
        end
    else
        allaics = result.cluster.aics;
        allbics = result.cluster.bics;
    end
    if isfield(result,'cluster')
        aics(end+1) = min(allaics);
        bics(end+1) = min(allbics);
    else
        warning([resultdir filesep files{i} ...
            ' is not a valid clustering file']);
    end
%     nclust(i) = ml_getfilenum(files{i});
    nclust(end+1) = result.cluster.bestk;
end
