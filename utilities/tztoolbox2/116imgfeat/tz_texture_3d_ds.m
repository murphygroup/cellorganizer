function feats = tz_texture_3d_ds(res,maskpath,ext,loadthre,medfilter)
%TZ_TEXTURE_3D_DS Obsolete.

%function feats = tz_texture_3d_ds(res,maskpath,ext,loadthre,medfilter)
%
%OVERVIEW:
%   calculate testure features for images under ds?
%PARAMETERS:
%   res - resource 

%This function is wrong! 22-MAY-2004 TINGZ

%Last modified by tingz on Mar. 26, 2003
%Last modified by tingz on Apr. 1, 2003     
%Last modified by tingz on Dec. 13, 2003

%error('This function is wrong. Please wait for repairing');

error(tz_genmsg('of','tz_texture_3d_ds'));

%load 3D image
image=xc_loadimage(res,ext,loadthre);

no_of_slices = size( image, 3);

%denoise
if medfilter == 1
    'median'
    for slice_no = 1:no_of_slices
        slice = image(:,:,slice_no);
        %denoise
        image(:,:,slice_no)=medfilt2(slice,[3,3]);
    end
end

%Change the type into double for calculation
image=double(image);

%size of the image
size(image)

%Get the max value of the image
max_pixel = max(image(:));
max_pixel
%return

%Scale the image into 0-255
image=image/(max_pixel/255);

%Change the type into uint8
image=uint8(image);

%Sub the background
image = mv_3dbgsub(image);

if( exist(maskpath) )
     load ( maskpath );
     maskpath
%     no_of_slices = size( image, 3);
     for slice_no = 1:no_of_slices
        slice = image(:,:,slice_no);
        slice(find(mask==0))=0;
        
        %s1=size(find(image(:,:,slice_no)==0));
        image(:,:,slice_no) = slice;
        %s2=size(find(image(:,:,slice_no)==0));
    end
end 

%Get the threshold of the image to binarize it
protthresh = 255*mb_nihthreshold( image);

%Binarize
protbin = mv_binarize( image, uint8(protthresh));

protimg=xc_mask(image,protbin);

feats=zeros(1,26);
z = xc_3Dtexture(protimg);
for m= 1:13
    feats(m) = z(m,14);
    feats(13+m) = z(m,15);
end

feats
