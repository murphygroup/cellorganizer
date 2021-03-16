function [features,names]=tz_imgwavobjfeat_mcf(rootdir,nobj,wavname)
%TZ_IMGWAVOBJFEAT_MCF Unknown.

%function [feats,names]=tz_imgwavobjfeat_mcf(rootdir,nobj,wavname)
%
%OVERVIEW:
%   calculate object features for mcf
%PARAMETERS:
%   rootdir - root directory
%   nobj - number of objects in each image
%   wavname - wavelet basis
%RETURN:
%   features - mcf
%   names - feature names
%DESCRIPTION:
%   
%HISTORY:
%   11-NOV-2004 Initial write TINGZ

classes=tz_cleandirs(mv_dir(rootdir));

for i=1:length(classes)
    imgdir=[rootdir '/' classes{i}];
    protdir=[imgdir '/prot'];
    dnadir=[imgdir '/dna'];
    cropdir=[imgdir '/crop'];
    protfiles=tz_cleandirs(mv_dir(protdir));
    dnafiles=tz_cleandirs(mv_dir(dnadir));
    cropfiles=tz_cleandirs(mv_dir(cropdir));
    cellfeats={};
    for j=1:length(protfiles)
        protimg=mv_readimage([protdir '/' protfiles{j}]);
        dnaimg=mv_readimage([dnadir '/' dnafiles{j}]);
        cropimg=mv_readimage([cropdir '/' cropfiles{j}]);
        objects=tz_wavfindobj(protimg,cropimg,nobj,wavname);
        dnaproc=ml_preprocess( dnaimg, cropimg, 'ml');
        feat=[];
        for k=1:nobj
            [feat(k,:),names]=tz_objfeat(objects{k},dnaproc);
        end
        cellfeats{j}=feat;
    end
    features{i}=cellfeats;
end

