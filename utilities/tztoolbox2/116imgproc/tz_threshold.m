function t = tz_threshold(img)
%TZ_THRESHOLD Thresholding uint8 image
%   T = TZ_THRESHOLD(IMG)
%   
%   See also

%   26-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

if strcmp(class(img),'uint8')~=1
    error('Wrong data type! The image should be uint8.');
end

MAXVALUE = 255;
img(img==0) = [];
img2 = uint8(MAXVALUE-double(img));
t1 = ml_rcthreshold(img2);
if sum(img2<t1) < 10
    img2(img2<t1)=[];
    t1 = ml_rcthreshold(img2);
end

t = MAXVALUE-t1+1;