function img2=tz_reorgimg(img,option)
%TZ_REORGIMG Rearrange the border of an image.
%   IMG2 = TZ_REORGIMG(IMG,OPTION)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2=tz_reorgimg(img,option)
%   
%OVERVIEW:
%   redefine the border of the image
%PARAMETERS:
%   img - input image
%RETURN:
%   img2 - output image
%
%DESCRIPTION:
%  H --> +
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

if ~exist('option','var')
    option='grv';
end

sc1=tz_distrforce(img,option);
sc2=tz_distrforce(img',option);

[minsc1,i1]=min(sc1);
[minsc2,i2]=min(sc2);

img2=img;

if i1<size(img,1)
    img2=[img((i1+1):end,:);img(1:i1,:)];
end

if i2<size(img,2)
    img2=[img2(:,(i2+1):end),img2(:,1:i2)];
end