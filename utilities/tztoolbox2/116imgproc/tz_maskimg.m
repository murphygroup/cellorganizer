function mask=tz_maskimg(img,imgtitle)
%TZ_MASKIMG Manually crop an image.
%   MASK = TZ_MASKIMG(IMG,IMGTITLE)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function mask=tz_maskimg(img)
%
%OVERVIEW:
%   build a mask for an image
%PARAMETERS:
%   img - input image
%   imgtitle - title when drawing mask
%RETURN:
%   mask - output mask
%DESCRIPTION:
%   mask either 2d image or 3d image
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments

if ~exist('imgtitle','var')
    imgtitle=[];
end

if size(img,3)>1
    img=max(img,[],3);
end

image_size = size( img);
imagesc( img);
truesize( image_size/1.5);
title(imgtitle);
mask = roipoly;