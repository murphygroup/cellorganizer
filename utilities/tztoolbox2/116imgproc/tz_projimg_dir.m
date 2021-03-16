function proj=tz_projimg_dir(imgdir,ext,maskfile,option)
%TZ_PROJIMG_DIR Project images under a directory.
%   TZ_PROJIMG_DIR(IMGDIR,EXT,MASKFILE,OPTION)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function proj=tz_projimg_dir(imgdir,ext,maskfile,option)
%
%OVERVIEW:
%   Get a projection of the images under a directioy
%PARAMETERS:
%   imgdir - the directory where the images are located
%   ext - extension of the image files. Any other file with a different extention will be ignored.
%   maskfile - a binary image file for masking
%   option - the way to get a projected image
%       'max'   pick the brightest pixel along z axis
%       'mean'  take mean of pixels along z axis
%       'sum'   take sum of pixels along z axis
%RETURN:
%   proj - the 2d image from projection
%DESCRIPTION:
%   There is no big difference between option 'mean' and 'sum'
%HISTORY:
%   ??-??-???? Initial write TINGZ
%   18-AUG-2004 Modified TINGZ
%       - add comments
%       - deal with empty file name
%
%See also tz_projimg_3d, tz_procimages_ds, tz_projimg_mcf.

if ~isempty(maskfile)
    mask=imread(maskfile);
else
    mask=[];
end

img=tz_loadimage(imgdir,ext,[]);

proj=tz_projimg_3d(img,mask,option);