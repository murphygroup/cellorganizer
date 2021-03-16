function img2 = tz_imcropbg(img,threshold)
%TZ_IMCROPBG Crop part of background out (Obsolete)
%   IMG2 = TZ_IMCROPBG(IMG) returns an image by cropping background of IMG.
%   Here background means pixels with value no greater than 0.
%   
%   IMG2 = TZ_IMCROPBG(IMG,THRESHOLD) defines background by a specified
%   value THRESHOLD.
%   
%   See also TZ_IMCROPBG2

%   20-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_imcropbg','ml_imcropbg'));

if nargin < 1
    error('1 or 2 arguments are required')
end

if nargin < 2
    threshold = 0;
end

img2=img;
img(img>threshold) = 1;
img(img<=threshold) = 0;

img2(:,sum(img,1)==0)=[];
img2(sum(img,2)==0,:)=[];