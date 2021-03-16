function tz_objdistimg_mcf(rootdir,savedir)
%TZ_OBJDISTIMG_MCF Distance maps for MCF images.
%   TZ_OBJDISTIMG_MCF(ROOTDIR,SAVEDIR)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_objdistimg_mcf(rootdir,savedir)
%OVERVIEW
%   generate distance maps for all images
%PARAMETERS
%   rootdir - image dir
%   savedir - saving dir
%RETURN
%
%DESCRIPTION
%   
%HISTORY
%   12-Apr-2005 Initial write TINGZ
%SEE ALSO
%   
class_names=tz_cleandirs(mv_dir(rootdir));
imgsize=[1024,1024];

if ~exist(savedir,'dir')
    unix(['mkdir ' savedir]);
end

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;

for class = first_class : last_class
    class_name = class_names{class};
    status=mkdir(savedir,class_name);
    saveclassdir=[savedir '/' class_name];
    
    status=mkdir(saveclassdir,'celldistmap');
    status=mkdir(saveclassdir,'dnadistmap');
    
    savecelldir=[saveclassdir '/celldistmap'];
    savednadir=[saveclassdir '/dnadistmap'];
    
    celldir=[rootdir '/' class_name '/cellbody'];
    dnadir=[rootdir '/' class_name '/nucbody'];
    
    cells=tz_cleandirs(mv_dir(celldir));
    ncell=length(cells);
    for N=1:ncell
        if ~exist([savecelldir '/' cells{N}],'file')
            load([celldir '/' cells{N}]);
            distimg=tz_objdistimg(cellbody,imgsize);
            save([savecelldir '/' cells{N}],'distimg');
        end
    end
    
    dnas=tz_cleandirs(mv_dir(dnadir));
    ndna=length(dnas);
    for N=1:ndna
        if ~exist([savednadir '/' dnas{N}],'file')
            load([dnadir '/' dnas{N}]);
            distimg=tz_objdistimg(nucbody,imgsize);
            save([savednadir '/' dnas{N}],'distimg');
        end
    end
end