function tz_kmeansaicbic_p(featfile,savedir,nclust)
%TZ_KMEANSAICBIC_P Cluster features.
%   TZ_KMEANSAICBIC_P(FEATFILE,SAVEDIR,NCLUST) cluster features into NCLUST 
%   clusters based on the features saved in FEATFILE. 10 trials will be
%   tried and aics and bics will be calculated. The variable in FEATFILE is
%   a [feature matrix] with the name features. SAVEDIR is the directory for
%   saving results.
%   
%   See also

%   05-Mar-2006 Initial write T. Zhao
%   02-May-2006 Modified T. Zhao
%       - ignore existing result file
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if license('checkout','statistics_toolbox')==0
    return
end

nCluster = nclust;
resultDirectory = savedir;
resultFile = ['cluster' num2str(nCluster)];
resultPath = [resultDirectory filesep resultFile '.mat'];
if exist(resultPath,'file')
    return;
end

addpath('/home/tingz/matlab/shared');
tz_initpath

data = load(featfile);

out = tz_clustering(data.features,struct('method','kmeansaic','ks', ...
    nclust,'tr',10,'zscore',1));

comments = ['clustering features in ' featfile];
script = ['tz_kmeansaicbic_p(''' ...
    featfile ''',''' savedir ''',' num2str(nclust) ')'];
version = 'debug';

tz_save(resultPath,{out,version,featfile},...
    {'cluster','version'},script,comments);