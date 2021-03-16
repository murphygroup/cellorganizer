function objects=tz_finddnaobjs_mcf(dirname,way)
%TZ_FINDDNAOBJS_MCF Obsolete.

%function objects=tz_finddnaobjs_mcf(dirname,way)
%
%OVERVIEW:
%   find objects in dna images under a directory
%PARAMETERS:
%   dirname - root directory
%   way - 'ml' or 'mb' for preprocessing
%RETURN:
%   objects - object cell array
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - revise comments

load('/home/velliste/feat/class_names.mat');

NumberOfClasses = 10;
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = mv_classname( class);
    protfiles = mv_dir([dirname '/' class_name '/dna/*.dat']);
    cropfiles = mv_dir([dirname '/' class_name '/crop/*.tif*']);
    L = length( protfiles);
    if( (L ~= length( cropfiles)))
        error( 'There must be a matching number of protein images, DNA images and crop  images'); 
    end
    n_images = L;
    cellobjs={};
    for N = 1 : n_images
        fprintf(1,[ class_name ', image ' num2str(N) '\n']);
        imagename = [dirname '/' class_name '/dna/' protfiles{N}];
        cropimagename = [dirname '/' class_name '/crop/' cropfiles{N}];
        
        cellobjs{N}=tz_findimgobjs(imagename,cropimagename,way);
    end
    objects{class}=cellobjs;
end