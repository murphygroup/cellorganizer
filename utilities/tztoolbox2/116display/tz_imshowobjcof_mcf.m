function tz_imshowobjcof_mcf(dirname,combobjects,combclass,combcellidx)
%TZ_IMSHOWOBJCOF_MDF Under construction

%function tz_imshowobjcof_mcf(dirname,combobjects,combclass,combcellidx)
%
%OVERVIEW:
%   show objects cof
%PARAMETERS:
%   dirname - directory for images
%   combobjects - combined objects
%   combclass - combined class label
%   combcellidx - combine cell idx
%RETURN:
%
%DESCRIPTION:
%   This funcion need updating
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%       - change function name tz_imshowdir_mcf -->tz_imshowobjcof_mcf

load('/home/velliste/feat/class_names.mat');

NumberOfClasses = 10;
first_class = 1;
last_class = NumberOfClasses;
for class = 2%first_class : last_class
    class_name = mv_classname( class);
    protfiles = mv_dir([dirname '/' class_name '/prot/*.dat']);
    cropfiles = mv_dir([dirname '/' class_name '/crop/*.tif*']);
    dnafiles = mv_dir([dirname '/' class_name '/dna/*.dat']);
    L = length( protfiles);
    if( (L ~= length( cropfiles)))
        error( 'There must be a matching number of protein images, DNA images and crop  images'); 
    end
    n_images = L;
    for N = 23%1 : n_images
        fprintf(1,[ class_name ', image ' num2str(N) '\n']);
        protimgname = [dirname '/' class_name '/prot/' protfiles{N}];
        cropimgname = [dirname '/' class_name '/crop/' cropfiles{N}];
        dnaimgname = [dirname '/' class_name '/dna/' dnafiles{N}];
        
        protimg=mv_readimage(protimgname);
        cropimg=mv_readimage(cropimgname);
        dnaimg=mv_readimage(dnaimgname);

        imshow(dnaimg,[])
        title(protimgname);
        hold on
        tz_plotgraph(tz_getobjcenter(combobjects(combclass==class & combcellidx==N)),[],1:2,'x','.',[])
        hold off
        drawnow
        pause(1)
        %cellobjs{N}=tz_findimgobjs(imagename,cropimagename,way);
    end
    %objects{class}=cellobjs;
end

% hold off
