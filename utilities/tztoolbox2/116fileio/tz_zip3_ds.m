function tz_zip3_ds( dirname1, task)
%TZ_ZIP3_DS Compress or uncompress images.
%   TZ_ZIP3_DS(DIRNAME,TASK) processes ds directory DIRNAME according
%   to the option of TASK:
%       'remove_nonimg' - remove non-tif images
%       'unzip' - uncompress images
%       'zip' - compress images
%   See also

%   ??-???-???? Initial write TINGZ
%   22-MAY-2004 Modified TINGZ
%       - remove return value
%       - change function name
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

prot_list = ml_dir(dirname1);
prot_list = prot_list(3:end);
for prot_no = 1:length(prot_list)
     prot_name = prot_list{prot_no};
     prot_fullpath = [dirname1 '/' prot_name];
     drug_list = ml_dir(prot_fullpath);
     drug_list = drug_list(3:end);
     ncondit = length(drug_list);
     ndrugs = ncondit/2;
     for drug_no = 1:ncondit
         drug_name = drug_list{drug_no};
         drug_fullpath = [prot_fullpath '/' drug_name];
         process_zip( drug_fullpath, task);
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function process_zip( dirname1, task)

image_names = ml_dir(dirname1);
image_names = image_names(3:end);
no_of_images = length(image_names);
for image_no = 1:no_of_images
    image_name = image_names{image_no};
    image_fullname = [dirname1 '/' image_name];
    slice_names = ml_dir(image_fullname);
    slice_names = slice_names(3:end);
    for s = 1:length(slice_names)
     slice_name = slice_names{s};
     switch( task)
      case 'remove_nonimg',
         if( length( findstr( slice_name, '.tif.')) > 0 | slice_name(end-3:end) == '.tif')
         else
             cmd = ['rm ' image_fullname '/"' slice_name '"']
             unix( cmd);
         end
      case 'unzip'
         if( slice_name(end-3:end) == '.bz2')
             cmd = ['bunzip2 ' image_fullname '/' slice_name]
             unix( cmd);
         end
     case 'zip'
         if( length( findstr( slice_name, '.tif.')) > 0 | slice_name(end-3:end) == '.tif')
             cmd = ['bzip2 ' image_fullname '/' slice_name]
             unix( cmd);
         end
      otherwise,
         error('Incorrect task name!');
     end
    end
end
