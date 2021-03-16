function img_out =tz_maskimg_3d(img,maskimg)
%TZ_MASKIMG_3D Crop 3D image by a 2D or 3D mask.
%   IMG2 = TZ_MASKIMG_3D(IMG,MASKIMG)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
%   10/1/14 Rewritten by Gregory Johnson to account for 2D masks on 3D
%   images, as well as made more efficient

%function img2=tz_maskimg_3d(img,maskimg)
%
%OVERVIEW:
%   

maskimg = maskimg > 0;

if size(img,3) ~= size(maskimg,3) & size(maskimg,3) == 1
    maskimg = repmat(maskimg, [1,1,size(img,3)]);
end

img_out = img.*maskimg;

