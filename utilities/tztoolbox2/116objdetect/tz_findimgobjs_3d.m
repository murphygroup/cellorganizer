function objects=tz_findimgobjs_3d(protimgdir,cropimgfile,loadthre)
%TZ_FINDIMGOBJS_3D Obsolete. See ML_FINDIMGOBJS_3D

%function objects=tz_findimgobjs_3d(protimgfile,cropimgfile,loadthre)
%   
%OVERVIEW:
%   find objects in 3d image files
%PARAMETERS:
%   protimgdir - image directory
%   cropimgfile -  mask image
%   loadthre - loading threshold
%RETURN:
%   objects - object cell array
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments

img = tz_loadimage( protimgdir, 'tif', loadthre);

img=double(img);
max_pixel = max(img(:));
img=img/(max_pixel/256);
img=uint8(img);
img = mv_3dbgsub(img);

mask=imread(cropimgfile);

% Mask the images
%    mask_img_fullpath = [dirname2 '/' image_name '.mat']
%    load ( mask_img_fullpath);

no_of_slices = size( img, 3);
for slice_no = 1:no_of_slices
    slice = img(:,:,slice_no);
    slice(find(mask==0))=0;
    img(:,:,slice_no) = slice;
end

protclean = img;

protthresh = 255*mb_nihthreshold( protclean);
protbin = mv_binarize( protclean, uint8(protthresh));
objects=mv_3dfindobj_sa( protbin,1);

nobj=length(objects);

for i=1:nobj
    nvoxel=size(objects{i}.voxels,2);
    objects{i}.voxels=objects{i}.voxels';
    sobj=[];
    idx=tz_index2order(objects{i}.voxels,size(protclean));
    objects{i}.voxels=[objects{i}.voxels,double(protclean(idx))];
end


