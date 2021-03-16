function tz_projimg_mcf(dirname,savedir,ext,option)
%TZ_PROJIMG_MCF Project MCF images.
%   TZ_PROJIMG_MCF(DIRNAME,SAVEDIR,EXT,OPTION)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_projimg_mcf(dirname,savedir,ext,option)
%
%OVERVIEW:
%   project mcf images
%PARAMETERS:
%   dirname - rootdir
%   savedir - saving directory
%   ext - extension of the image files
%   option - projection option. see tz_projimg_dir
%RETURN:
%
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   02-NOV-2004 Modified TINGZ
%       - add comments

class_names=tz_cleandirs(mv_dir(dirname));

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
    
    cells=tz_cleandirs(mv_dir([dirname '/' class_name]));
    ncell=length(cells);
    for N=1:ncell
       % status=mkdir(saveclassdir,cells{N});
        %savecelldir=[saveclassdir '/' cells(N)];
        celldir=[dirname '/' class_name '/' cells{N}];
        protdir = [celldir '/prot']
        cropfiles = mv_dir([celldir '/crop/*.tif*'])
        proj=tz_projimg_dir(protdir,ext,[celldir '/crop/' cropfiles{1}],option);
        save([saveclassdir '/' cells{N} '.mat'],'proj')
    end
end

