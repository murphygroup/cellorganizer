function dimg = tz_imsub(img1,img2,offset,mask)
%TZ_IMSUB Subtraction between two images.
%   DIMG = TZ_IMSUB(IMG1,IMG2) returns IMG1-IMG2. DIMG is a matrix of
%   double.
%
%   DIMG = TZ_IMSUB(IMG1,IMG2,OFFSET) takes the offset of IMG2. This means
%   the original point [0 0] in IMG2 has OFFSET in IMG1.
%
%   DIMG = TZ_IMSUB(IMG1,IMG2,OFFSET,MASK) takes the mask of IMG2.
%
%   See also

%   24-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('At least 2 arguments are required')
end

img2 = tz_impad2(img2,img1);

if nargin >= 3
    img2 = tz_imtranslate(img2,offset);
end   
 
if nargin>=4
    mask = tz_impad2(mask,img1);
    mask = tz_imtranslate(mask,offset);
end

dimg = double(img1)-double(img2(1:size(img1,1),1:size(img1,2)));

if exist('mask','var')
    dimg(mask==0) = 0;
end