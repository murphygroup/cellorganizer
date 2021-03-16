function objects=tz_findmatobjs_mcf(dirname,way)

%function objects=tz_findimgobjs_mcf
%
%OVERVIEW:
%   find objects in images under a directory
%PARAMETERS:
%   dirname - directory
%   way - 'ml' or 'mb' for preprocessing
%RETURN:
%   objects - object cell array
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments

mv_classname=tz_cleandirs(mv_dir(dirname));

NumberOfClasses = length(mv_classname);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = mv_classname{class};
    protfiles = mv_dir([dirname '/' class_name '/*.mat']);
    %cropfiles = mv_dir([dirname '/' class_name '/crop/*.tif*']);
    L = length( protfiles);
%     if( (L ~= length( cropfiles)))
%         error( 'There must be a matching number of protein images, DNA images and crop  images'); 
%     end
    n_images = L;
    cellobjs={};
    for N = 1 : n_images
        fprintf(1,[ class_name ', image ' num2str(N) '\n']);
        imagename = [dirname '/' class_name '/' protfiles{N}];
%         cropimagename = [dirname '/' class_name '/crop/' cropfiles{N}];
        
        cellobjs{N}=tz_findmatobjs(imagename,way);
    end
    objects{class}=cellobjs;
end