function tz_maskdir_mcf(dirname,ext)
%TZ_MASKDIR_MCF Crop MCF images.
%   TZ_MASKDIR_MCF(DIRNAME,EXT)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_maskdir_mcf(dirname,ext)
%
%OVERVIEW:
%
%PARAMETERS:
%   dirname - root directory
%   ext - extension of image files
%RETURN:
%
%DESCRIPTION:
%   only for 3d images
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments
%

class_names=tz_cleandirs(mv_dir(dirname));

savedir=dirname;

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = class_names{class};
    saveclassdir=[savedir '/' class_name];
    cells=tz_cleandirs(mv_dir([dirname '/' class_name]));
    ncell=length(cells);
    for N=1:ncell
        celldir=[dirname '/' class_name '/' cells{N}];
        protdir = [celldir '/cell'];
        cropdir = [celldir '/crop'];
        tz_maskdir(protdir,ext,cropdir,'3d');
    end
end
