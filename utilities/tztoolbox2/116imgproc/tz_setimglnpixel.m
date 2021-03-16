function img2=tz_setimglnpixel(img,s,t,val)
%TZ_SETIMGLNPIXEL Draw a line in an image.
%   IMG2 = TZ_SETIMGLNPIXEL(IMG,S,T) set all pixels on a line to 1 in the
%   image IMG. The line is from S to T. See ML_GETLINEPTS for S and T.
%   
%   IMG2 = TZ_SETIMGLNPIXEL(IMG,S,T,VAL) set the value VAL instead of 1.
%
%   See also TZ_SETIMGLNPIXEL2 ML_GETLINEPTS

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin<4
    val=1;
end

pts=ml_getlinepts(s,t);
img2=ml_setimgptspixel(img,[pts,zeros(size(pts,1),1)+val]);