function tz_kmeansaicbic2_p(featfile,paramfile,nclust)
%TZ_KMEANSAICBIC_P Cluster features.
%   TZ_KMEANSAICBIC_P(FEATFILE,PARAMFILE,NCLUST) cluster features into
%   based on the NCLUST features.
%   features saved in FEATFILE. 10 trials will be tried and aics and bics
%   will be calculated. The variable in FEATFILE is a [feature matrix] with
%   the name features. SAVEDIR is the directory for saving results.
%   
%   See also TZ_KMEANSAICBIC_P

%   29-Mar-2006 Initial write T. Zhao
%   01-JUN-2006 Modified T. Zhao
%      - fix the bug of number of clusters
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if license('checkout','statistics_toolbox')==0
    return
end

param = load(paramfile);

nCluster = nclust;
resultDirectory = param.savedir;
resultFile = ['cluster' num2str(nCluster)];
resultPath = [resultDirectory filesep resultFile '.mat'];
if exist(resultPath,'file')
    return;
end

addpath('/home/tingz/matlab/shared');
tz_initpath

data = load(featfile);

param.clustering.ks = nclust;
out = tz_clustering(data.features,param.clustering);

comments = ['clustering features in ' featfile];
script = ['tz_kmeansaicbic2_p(''' ...
    featfile ''',''' paramfile ''',''' num2str(nclust) ''')'];
version = 'debug';

tz_save(resultPath,{out,version,featfile},...
    {'cluster','version'},script,comments);
