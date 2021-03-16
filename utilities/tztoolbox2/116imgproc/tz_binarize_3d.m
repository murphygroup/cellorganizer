function tz_binarize_3d(res,des,maskpath,ext,loadthre,imgtype,med)
%TZ_BINARIZE_3D Binarize 3d images.
%   TZ_BINARIZE_3D(RES,DES,MASKPATH,EXT,LOADTHRE,IMGTYPE,MED)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_binarize_3d(res,des,maskpath,ext,loadthre,imgtype,med)
%OVERVIEW:
%   Binarize 3d images with extension ext and save the results as imgtype under des
%PARAMETERS:
%   res - resouce directory
%   des - distination directory
%   maskpath - full path for the mask
%   ext - resource file extension
%   loadthre - threshold to remove noise
%   imgtype - saving file type
%   med - median filtering or not
%RETURN:
%   no return
%DESCRIPTION:
%   The thresholding is automatic
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   06-MAR-2003 Modified TINGZ
%   12-MAR-2003 Modified TINGZ
%     -- Change path
%   01-APR-2003 Modified TINGZ
%   07-JUL-2003 Modified TINGZ
%   12-JUL-2003 Modified TINGZ
%   16-JUL-2003 Modified TINGZ
%   13-DEC-2003 Modified TINGZ

image = xc_loadimage( res, ext, loadthre);

res

no_of_slices = size( image, 3);

if med==1
%denoise
'medfilter'
for slice_no = 1:no_of_slices
    slice = image(:,:,slice_no);
    %denoise
    image(:,:,slice_no)=medfilt2(slice,[3,3]);
end
end

% Convert image to double and substract the background
image=double(image);
max_pixel = max(image(:));
image=image/(max_pixel/256);
image=uint8(image);
image = ml_3dbgsub(image);

if( exist(maskpath) )
     load ( maskpath );
     maskpath

     for slice_no = 1:no_of_slices
        slice = image(:,:,slice_no);
        slice(find(mask==0))=0;
        
        image(:,:,slice_no) = slice;
     end
end 

%Get the threshold of the image to binarize it
protthresh = 255*mb_nihthreshold( image);

%Binarize
protbin = mv_binarize( image, uint8(protthresh));

%Change the result into a gray image
grayimg=mat2gray(double(protbin));

%maximg=max(grayimg(:))

%Save
i='0';
for slice_no = 1:no_of_slices
    desfilename=['0' i+floor(slice_no/10) i+slice_no-floor(slice_no/10)*10 '.' imgtype];
    [des '/' desfilename]
    
    imwrite(grayimg(:,:,slice_no),[des '/' desfilename],imgtype);
end
