function tz_makemasks(diroot,dpdir,dmdir)
%TZ_MAKEMASKS Make masks.
%   TZ_MAKEMASKS(DIROOT,DPDIR,DMDIR)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_makemasks(diroot,dpdir,dmdir)
%
%OVERVIEW:
%   make masks
%PARAMETERS:
%   diroot - root directory
%   dpdir - protein directory
%   dmdir - mask diretory
%RETURN:
%
%DESCRIPTION:
%   This function is isolated from process_all to make masks in the computer 
%   without original images.
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

prot_list = ml_dir([diroot '/' dpdir])
prot_list = prot_list(3:end)

for prot_no = 1:length(prot_list)
     prot_name = prot_list{prot_no}
     prot_fullpath = [diroot '/' dpdir '/' prot_name];
     mask_prot_fullpath = [diroot '/' dmdir '/' prot_name];
     
     make_dir( mask_prot_fullpath);
     
     drug_list = ml_dir(prot_fullpath);
     drug_list = drug_list(3:end);
     ncondit = length(drug_list);
     ndrugs = ncondit/2;
     
     for drug_no = 1:ncondit
         drug_name = drug_list{drug_no};
         drug_fullpath = [prot_fullpath '/' drug_name];
        
         mask_fullpath = [mask_prot_fullpath '/' drug_name];
         make_dir( mask_fullpath);
         make_mask(drug_fullpath,mask_fullpath);
               
     end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_mask( dirname1, dirname2)

image_names = ml_dir(dirname1);
image_names = image_names(3:end);

no_of_images = length(image_names);
feature_matrix = [];
for image_no = 1:no_of_images
    image_name = image_names{image_no};
    image_fullname = [dirname1 '/' image_name]
    mask_img_fullpath = [dirname2 '/' image_name];
    
    if( ~exist( mask_img_fullpath))
        % Read projection image
        load( image_fullname);
        image_size = size( proj_image);
        imagesc( proj_image);
        truesize( image_size/1.5);
        mask = roipoly;
        save( mask_img_fullpath, 'mask') ;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function make_dir( dirname)

if( ~exist( dirname, 'dir'))
    command = ['mkdir ' dirname]
    status = unix( command);
    if( status ~= 0)
        error(['Cannot create directory: ' dirname]);
    end
end
