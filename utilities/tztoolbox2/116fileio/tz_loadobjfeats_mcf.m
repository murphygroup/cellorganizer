function features = tz_loadobjfeats_mcf(featdir)
%TZ_LOADFEATS_MCF Load object-level features.
%   FEATURES = TZ_LOADFEATS_MCF(FEATDIR) loads object-level features from
%   directory FEATDIR into a cell array FEATURES. Usually FEATDIR is a
%   directory containing subdirectories, which represent proteins. Feaure
%   files are saved under each protein subdirectory. Each feature file
%   stores features of objects in one cell. So FEATURES{i}{j} is the
%   features of objects in the jth cell and ith class. Each row contains
%   features of one object.

%   29-NOV-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

classes=tz_cleandirs(mv_dir(featdir));

nclass=length(classes);

for i=1:nclass
    featfiles=tz_cleandirs(mv_dir([featdir '/' classes{i}]));
    ncell=length(featfiles);
   
    for j=1:length(featfiles);
        load([featdir '/' classes{i} '/' featfiles{j}])
        features{i}{j}=feats;
    end
end
