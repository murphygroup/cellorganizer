function [maics,ks]=tz_loadaics(dirname,task)

%TZ_LOADAICS Load AIC values from files.
%   MAICS = TZ_LOADAICS(DIRNAME,TASK) loads AIC values of clustering
%   from files under directory DIRNAME. TASK must be 'pl' currently. 
%
%   [MAICS,KS] = TZ_LOADAICS(DIRNAME,TASK) also returns KS, which
%   is a vector of numbers of clusters
%
%   See also

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

files=mv_dir(dirname);
files=files(3:end);
maics=[];
ks=[];
for i=1:length(files)
    fullpath=[dirname '/' files{i}]
    load(fullpath);
    maics=[maics;aics];
    ks(i)=tz_getfilenum(files{i});
    switch task
    case 'pl'
        for i=1:length(posts(:))
            posts{i}=tz_post2label(posts{i});
        end
        save(fullpath,'aics','centers','posts');
    end
end