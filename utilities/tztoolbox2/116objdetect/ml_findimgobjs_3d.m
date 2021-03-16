function objects = ml_findimgobjs_3d(protimgdir,cropimgfile,loadthre)
%ML_FINDIMGOBJS_3D Find objects in 3D image files.
%   OBJECTS = ML_FINDIMGOBJS_3D(PROTIMGDIR,CROPIMGFILE) returns a cell
%   array of 3D objects and holes which are extracted from the 3D image
%   files under the directory PROTIMGDIR. The files under PROTIMGDIR should
%   be the slices of the same 3D image.
%   
%   OBJECTS = ML_FINDIMGOBJS_3D(PROTIMGDIR,CROPIMGFILE,LOADTHRE) remove
%   bright pixels from the 3D image before finding objects, i.e. all values
%   greater than LOADTHRE well be set to 0.
%   
%   See also

%   25-Oct-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu

if ~exist('loadthre','var')
    loadthre = [];
end

%img = ml_loadimage( protimgdir, 'tif', loadthre);

%img=double(img);
img = protimgdir;
max_pixel = max(img(:));
img=img/(max_pixel/256);
img=uint8(img);
img = ml_3dbgsub(img);

% mask=imread(cropimgfile);
mask = cropimgfile;

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

protthresh = 255*ml_threshold( protclean);
protbin = ml_binarize( protclean, uint8(protthresh));
objects=ml_3dfindobj_sa( protbin,1);

nobj=length(objects);

for i=1:nobj
    nvoxel=size(objects{i}.voxels,2);
    objects{i}.voxels=objects{i}.voxels';
    sobj=[];
    %idx=tz_index2order(objects{i}.voxels,size(protclean));
    idx = sub2ind(size(protclean),objects{i}.voxels(:,1), ...
            objects{i}.voxels(:,2),objects{i}.voxels(:,3));
    objects{i}.voxels=[objects{i}.voxels,double(protclean(idx))];
end


