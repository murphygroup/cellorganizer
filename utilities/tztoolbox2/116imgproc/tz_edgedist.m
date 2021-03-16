function img2=tz_edgedist(img)
%TZ_EDGEDIST Build a distance map for the edge of an image.
%   IMG2 = TZ_EDGEDIST(IMG) returns a distance map of the edge in IMG, of
%   which pixels with intensity 0 are taken as background.
%   
%   See also

%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

img2=img>0;
img2=bwperim(img2);
img2=bwdist(img2);
