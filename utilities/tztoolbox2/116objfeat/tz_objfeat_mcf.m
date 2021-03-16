function [features,names]=tz_objfeat_mcf(objects,rootdir)
%TZ_OBJFEAT_MCF Under construction.

%function [feats,names]=tz_imgwavobjfeat_mcf(rootdir,nobj,wavname)
%
%OVERVIEW:
%   calculate object features for mcf
%PARAMETERS:
%   objects - mcf objects
%   rootdir - root directory for dna images
%RETURN:
%   features - mcf
%   names - feature names
%DESCRIPTION:
%   
%HISTORY:
%   11-NOV-2004 Initial write TINGZ

classes=tz_cleandirs(mv_dir(rootdir));

for i=1:length(objects)
    imgdir=[rootdir '/' classes{i}];
    
    dnadir=[imgdir '/dna'];
    cropdir=[imgdir '/crop'];
    
    dnafiles=tz_cleandirs(mv_dir(dnadir));
    cropfiles=tz_cleandirs(mv_dir(cropdir));
    cellfeats={};
    for j=1:length(objects{i})
        [i j]
        dnaimg=mv_readimage([dnadir '/' dnafiles{j}]);
        cropimg=mv_readimage([cropdir '/' cropfiles{j}]);
        
        dnaproc=ml_preprocess( dnaimg, cropimg, 'ml');
        feat=[];
        for k=1:length(objects{i}{j})
            [feat(k,:),names]=tz_objfeat(objects{i}{j}{k},dnaproc);
        end
        cellfeats{j}=feat;
    end
    features{i}=cellfeats;
end