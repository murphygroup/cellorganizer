function tz_imedge_mcf(rootdir,savedir)
%TZ_IMEDGE_MCF Extract and process edges for MCF cell and DNA channels.
%   TZ_IMEDGE_MCF(ROOTDIR,SAVEDIR)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_imedge_mcf(rootdir,savedir)
%OVERVIEW:
%   Batch processing of extracting edges for cells and nucleus
%PARAMETERS:
%   rootdir - root directory
%   savedir - save directory
%RETURN:
%   
%DESCRIPTION
%
%HISTORY:
%   27-SEP-2004 Initial write TINGZ
%

class_names=tz_cleandirs(mv_dir(rootdir));

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
    
    status=mkdir(saveclassdir,'celledge');
    status=mkdir(saveclassdir,'dnaedge');
    
    savecelldir=[saveclassdir '/celledge'];
    savednadir=[saveclassdir '/dnaedge'];
    
    celldir=[rootdir '/' class_name '/cell'];
    dnadir=[rootdir '/' class_name '/dna'];
    
    cells=tz_cleandirs(mv_dir(celldir));
    ncell=length(cells);
    for N=1:ncell
        if ~exist([savecelldir '/' cells{N}],'file')
            load([celldir '/' cells{N}]);
            celedge=tz_imedge(selimg,[],'ce');
            imshow(celedge,[]);
            drawnow
            save([savecelldir '/' cells{N}],'celedge');
        end
    end
    
    dnas=tz_cleandirs(mv_dir(dnadir));
    ndna=length(dnas);
    for N=1:ndna
        if ~exist([savednadir '/' dnas{N}],'file')
            load([dnadir '/' dnas{N}]);
            dnaedge=tz_imedge(selimg,[],'nu');
            imshow(dnaedge,[]);
            drawnow
            save([savednadir '/' dnas{N}],'dnaedge');
        end
    end
end