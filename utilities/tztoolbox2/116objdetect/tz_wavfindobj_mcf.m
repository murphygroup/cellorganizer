function objects=tz_wavfindobj_mcf(rootdir,savepath,nobj,wavname, ...
    level,featset,alpha,option,t,dnat)
%TZ_WAVFINDOBJ_MCF Wavelet-based object detection for MCF images.
%   OBJECTS = TZ_WAVFINDOBJ_MCF(ROOTDIR,SAVEPATH,NOBJ,WAVNAME,LEVEL,
%   FEATSET,ALPHA,OPTION,T,DNAT)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [feats,names]=tz_wavfindobj_mcf(rootdir,savepath,nobj,wavname,level,featset,alpha,option,t)
%
%OVERVIEW:
%   calculate object features for mcf
%PARAMETERS:
%   rootdir - root directory
%   savpath - path for save
%   nobj - number of objects in each image
%   wavname - wavelet basis
%   savepath - path for saving
%   option - how to preprocess
%RETURN:
%   features - mcf
%   names - feature names
%DESCRIPTION:
%   
%HISTORY:
%   11-NOV-2004 Initial write TINGZ

classes=tz_cleandirs(mv_dir(rootdir));

if exist(savepath,'file')
    eobjs=load(savepath);
else
    eobjs.objects={};
end

if isempty(eobjs.objects)
    startclass=1;
else
    startclass=length(eobjs.objects);
end

objects=eobjs.objects;

for i=startclass:length(classes)
    imgdir=[rootdir '/' classes{i}];
    protdir=[imgdir '/prot'];
    dnadir=[imgdir '/dna'];
    cropdir=[imgdir '/crop'];
    protfiles=tz_cleandirs(mv_dir(protdir));
    dnafiles=tz_cleandirs(mv_dir(dnadir));
    cropfiles=tz_cleandirs(mv_dir(cropdir));
    obj={};
    if i>length(eobjs.objects)
        startcell=1;
    else
        startcell=length(eobjs.objects{i})+1;
    end
    
    for j=startcell:length(protfiles)
        protimg=mv_readimage([protdir '/' protfiles{j}]);
        dnaimg=mv_readimage([dnadir '/' dnafiles{j}]);
        cropimg=mv_readimage([cropdir '/' cropfiles{j}]);
        dnaproc=ml_preprocess(dnaimg, cropimg, 'ml', 'yes');
        obj=tz_wavfindobj(protimg,cropimg,dnaproc,dnat,nobj,wavname,level,featset,alpha,option,t);

        if j==1
            objects{i}{1}=obj;
        else
            objects{i}={objects{i}{1:end},obj};
        end
        save(savepath,'objects','nobj','wavname','level','featset','alpha','option','t','dnat');
    end
end

