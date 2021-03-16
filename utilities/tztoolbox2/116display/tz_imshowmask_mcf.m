function tz_imshowmask_mcf(dirname,ext)
%TZ_IMSHOWMASK_MCF Show 3D MCF images.
%   TZ_IMSHOWMASK_MCF(DIRNAME,EXT) shows images with extension EXT under 
%   the MCF directory DIRNAME.

%   16-Sep-2005 Initial write T. Zhao
%   ??-???-2004 Initial write T. Zhao
%   01-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

class_names=tz_cleandirs(mv_dir(dirname));
figure

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = class_names{class};
    
    cells=tz_cleandirs(mv_dir([dirname '/' class_name]));
    ncell=length(cells);
    for N=1:ncell
       % status=mkdir(saveclassdir,cells{N});
        %savecelldir=[saveclassdir '/' cells(N)];
        celldir=[dirname '/' class_name '/' cells{N}];
        protdir = [celldir '/prot'];
        cropfiles = mv_dir([celldir '/crop/*.tif*']);
        mask=imread([celldir '/crop/' cropfiles{1}]);
        img=tz_loadimage(protdir,ext,[]);
        img=max(img,[],3);
        tz_imshowmask(img,mask);
        title(celldir);
        celldir
    end
end