function tz_scanimages(dirname)
%TZ_SCANIMAGES Check format of images in ds.
%   TZ_SCANIMAGES(DIRNAME) check if the file types are consistent under
%   ds DIRNAME. This function provides a template for batch processing
%   ds.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

prot_list = ml_dir(dirname);
prot_list = prot_list(3:end);

for prot_no = 1:length(prot_list)
     prot_name = prot_list{prot_no};
     prot_fullpath = [dirname '/' prot_name];
      
     
     drug_list = ml_dir(prot_fullpath);
     drug_list = drug_list(3:end);
     ncondit = length(drug_list);
     
     for drug_no = 1:ncondit
         drug_name = drug_list{drug_no};
         drug_fullpath = [prot_fullpath '/' drug_name];
         checkincon( drug_fullpath);
     end
 end
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkincon(dirname)

image_names = ml_dir(dirname);
image_names = image_names(3:end);
no_of_images = length(image_names);

for image_no = 1:no_of_images
    image_name = image_names{image_no};
    image_fullname = [dirname '/' image_name]
    
    dir1 = ml_dir([image_fullname '/*.bz2']);
    dir2 = ml_dir([image_fullname '/*.tif']);
    if (length(dir1)~=0)&(length(dir2)~=0)
        image_fullname
        a='warning!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        pause
    end
end
