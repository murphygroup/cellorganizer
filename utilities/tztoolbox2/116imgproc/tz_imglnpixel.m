function ps=tz_imglnpixel(img,s,t)
%TZ_IMGLNPIXEL Get image pixel values along a line.
%   PS = TZ_IMGLNPIXEL(IMG,S,T)
%   
%   See also IMPROFILE, TZ_IMGLNPIXEL2

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function pts=tz_imglnpixel(img,s,t)

pts=tz_getlinept(s,t);
ps=tz_imgptspixel(img,pts);