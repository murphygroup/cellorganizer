function tz_kmeansaicbic4_p(featfile,paramfile,savedir,nclust)
%TZ_KMEANSAICBIC4_P Cluster features.
%   TZ_KMEANSAICBIC4_P(FEATFILE,PARAMFILE,SAVEDIR,NCLUST) cluster features
%   into NCLUST clusters based on features saved in FEATFILE. 10 trials
%   will be tried and aics and bics will be calculated. The variable in
%   FEATFILE is a [feature matrix] with the name features. PARAMFILE
%   contains the clustering parameters. SAVEDIR is the directory for saving
%   results. 
%   
%   See also TZ_KMEANSAICBIC_P TZ_KMEANSAICBIC2_P

%   29-Mar-2006 Initial write T. Zhao
%   01-JUN-2006 Modified T. Zhao
%      - fix the bug of number of clusters
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 4
    error('Exactly 4 arguments are required')
end

if license('checkout','statistics_toolbox')==0
    return
end

param = load(paramfile);

nCluster = nclust;
resultDirectory = savedir;
if ~exist(resultDirectory,'dir')
    mkdir(resultDirectory);
end

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
script = ['tz_kmeansaicbic4_p(''' ...
    featfile ''',''' paramfile ''',''' savedir ''',''' num2str(nclust) ''')'];
version = 'debug';

tz_save(resultPath,{out,version,featfile},...
    {'cluster','version'},script,comments);
