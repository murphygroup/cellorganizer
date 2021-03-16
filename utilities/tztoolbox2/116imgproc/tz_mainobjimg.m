function img2=tz_mainobjimg(img)
%TZ_MAINOBJIMG Remove all objects except the biggest one in an image.
%   IMG2 = TZ_MAINOBJIMG(IMG) returns an image that only contains the
%   biggest object in IMG.
%   
%   See also

%   ??-???-???? Initial write T. Zhao
%   01-NOV-2004 Modified T. Zhao
%       - add comments
%       - extended to gray level image
%       - change function name tz_mainobjimg_bw --> tz_mainobjimg
%   Copyright (c) Murphy Lab, Carnegie Mellon University


limg=bwlabel(img>0);
objnum=max(limg(:));
        
lhist=[];
for i=1:objnum
    lhist(i)=sum(sum(limg==i));
end
[y,maxl]=max(lhist);

img(limg~=maxl)=0;

img2=img;
