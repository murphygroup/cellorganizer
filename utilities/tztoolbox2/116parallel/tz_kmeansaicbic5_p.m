function tz_kmeansaicbic5_p(featfile,paramfile,savedir,nclusts)
%TZ_KMEANSAICBIC5_P Cluster features.
%   TZ_KMEANSAICBIC5_P(FEATFILE,PARAMFILE,SAVEDIR,NCLUSTS) cluster features
%   into NCLUSTS clusters based on features saved in FEATFILE. The variable in
%   FEATFILE is a [feature matrix] with the name features. PARAMFILE
%   contains the clustering parameters. It contains a variable called 
%   'clustering', which the structure parameter of the function TZ_CLUSTERING.
%   The field 'method' in the structure should be either 'kmeansaic' or 
%   'kmeansbic'. The field 'ks' has no use. The number of clusters is
%   specified by NCLUSTS instead. SAVEDIR is the directory for saving results. 
%   
%   See also TZ_KMEANSAICBIC_P TZ_KMEANSAICBIC2_P TZ_KMEANSAICBIC4_P

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

addpath('/home/tingz/matlab/shared');
tz_initpath
    
param = load(paramfile);

for nclust=nclusts
    nCluster = nclust;
    resultDirectory = savedir;
    if ~exist(resultDirectory,'dir')
        mkdir(resultDirectory);
    end
    
    resultFile = ['cluster' num2str(nCluster)];
    resultPath = [resultDirectory filesep resultFile '.mat'];
    
    [s,msg] = mkdir([resultPath '.ctr']);
    if strfind(msg,'exist')
        continue;
    end
        
    if exist(resultPath,'file')
        rmdir([resultPath '.ctr']);
        continue;
    end
    


    data = load(featfile);
    
    param.clustering.ks = nclust;
    out = tz_clustering(data.features,param.clustering);
    
    comments = ['clustering features in ' featfile];
    script = ['tz_kmeansaicbic5_p(''' ...
              featfile ''',''' paramfile ''',''' savedir ''',''' ...
              num2str(nclust) ''')'];
    version = 'debug';

    tz_save(resultPath,{out,version,featfile},...
            {'cluster','version'},script,comments);
    rmdir([resultPath '.ctr']);
end
