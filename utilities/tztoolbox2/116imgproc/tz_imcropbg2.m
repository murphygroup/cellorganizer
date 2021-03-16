function img2 = tz_imcropbg2(img,mask)
%TZ_IMCROPBG2 Crop part of background in a mask out
%   IMG2 = TZ_IMCROPBG2(IMG,MASK) returns an image by cropping background
%   out. Here background means pixels with value no greater than 0 in MASK.
%   
%   See also TZ_IMCROPBG

%   10-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

img2=img;

img2(:,sum(mask,1)==0)=[];
img2(sum(mask,2)==0,:)=[];