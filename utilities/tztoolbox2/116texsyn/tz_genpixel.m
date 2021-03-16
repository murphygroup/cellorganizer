function rpixel=tz_genpixel(img,mask,wnd,qmask,df,maxmasked,ep)

%function rpixel=tz_genpixel(img,mask,wnd,qmask,df,maxmasked,ep)
%
%OVERVIEW:
%   generate a pixel in the center of the window
%PARAMETERS:
%   

[mwnd,mmask,mpos,mdist]=tz_getbestmatch(img,mask,wnd,qmask,df,maxmasked);
cpixels=tz_getknnpixels(img,mask,mwnd,qmask,df,maxmasked,ep);
rpixel=cpixels(randperm(length(cpixels)));
rpixel=rpixel(1);