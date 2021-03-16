function [names,features] = tz_matobjfeats_mcf(dirname,featsetname,savedir)
%TZ_MATOBJFEATS_MCF Calculate object features for matlab data files.
%   NAMES = TZ_MATOBJFEATS_MCF(DIRNAME,FEATSETNAME,SAVEDIR)
%   
%   [NAMES,FEATURES] = TZ_MATOBJFEATS_MCF(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [names,features] = tz_matobjfeats_mcf(dirname,featsetname)
%   
%OVERVIEW:
%   calculate object features for matlab data files
%PARAMETERS:
%   dirname - root directory
%   featsetname - feature set
%   savedir - save directory
%RETURN:
%   names - feature names
%   features - object features
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments

if exist([savedir '/' '3d2dobjfeats.mat'],'file')
    load([savedir '/' '3d2dobjfeats.mat']);
else
    features={};
end

mv_classname=tz_cleandirs(mv_dir(dirname));

NumberOfClasses = length(mv_classname);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = mv_classname{class};
    protfiles = mv_dir([dirname '/' class_name '/prot/*.mat']);
    dnafiles = mv_dir([dirname '/' class_name '/dna/*.mat']);
    
    L = length( protfiles);

    n_images = L;
    for N = 1 : n_images
        if class>length(features)
            len=0;
        else
            len=length(features{class});
        end
        fprintf(1,[ class_name ', image ' num2str(N) '\n']);
        if N>len
            
            imagename = [dirname '/' class_name '/prot/' protfiles{N}];
            dnaimagename = [dirname '/' class_name '/dna/' dnafiles{N}];
            load(imagename);
            protimg=double(selimg);
            load(dnaimagename);
            dnaimg=double(selimg);
            cropimg=ones(size(protimg));
            
            [names, feats] = mv_objfeatset( protimg, cropimg, dnaimg, featsetname, 0);
            features{class}{N}=feats;
        end
        save([savedir '/' '3d2dobjfeats.mat']);
    end
    
end