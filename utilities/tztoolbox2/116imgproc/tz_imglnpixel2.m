function ps=tz_imglnpixel2(img,s,a,len)
%TZ_IMGLNPIXEL Get image pixel values along a line.
%   PS = TZ_IMGLNPIXEL(IMG,S,T)
%   
%   See also IMPROFILE, TZ_IMGLNPIXEL

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function pts=tz_imglnpixel2(img,a,len)

if len==-1
    len=sqrt(sum(size(img).^2));
end

pts=tz_getlinept2(s,a,len);

ps=tz_imgptspixel(img,pts);