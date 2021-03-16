function img2 = tz_tophat(img,wndsize)
%TZ_TOPHAT Top hat algorithm.
%   IMG2 = TZ_TOPHAT(IMG,WNDSIZE)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2 = tz_tophat(img,wndsize)
%OVERVIEW:
%   top hat algorithm
%PARAMETERS:
%   img - input image
%   wndsize - window size
%RETURN:
%   img2 - output image
%DESCRIPTION:
%
%HISTORY:
%   27-JUL-2003 Initial write TINGZ

[height width]=size(img)
img2=img;

if width==1
    loopj=1
else
    loopj=wndsize(2)+1:width-wndsize(2);
end

if height==1
    loopi=1;
else
    loopi=wndsize(1)+1:height-wndsize(1);
end

for j=loopj
    for i=loopi
        maskwnd=img(i-wndsize(1):i+wndsize(1),j-wndsize(2):j+wndsize(2));
        img2(i,j)=min(maskwnd(:));
    end
end

img=img2;

for j=loopj
    for i=loopi
        maskwnd=img(i-wndsize(1):i+wndsize(1),j-wndsize(2):j+wndsize(2));
        img2(i,j)=max(maskwnd(:));
    end
end

img2=img-img2;
