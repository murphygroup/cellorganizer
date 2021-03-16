function img2=tz_setimgptspixel(img,pts)
%TZ_SETIMGPTSPIXEL Obsolete.

%function img2=tz_setimgptspixel(img,pts)
%OVERVIEW:
%   set values in images
%PARAMETERS:
%   img - input image
%   pts - positions (values) nx2(3), if it has 3 columns,
%       the last column is the values. if not, all set to 1
%RETURN:
%   img2 - output image

error('Function tz_setimgptspixel is out of date. Please use ml_setimgptspixel');

img2=img;
imgsize=size(img);
if size(pts,2)<3
    pts(:,3)=1;
end

pts(pts(:,1)<=0 | pts(:,1)>imgsize(1),:)=[];
pts(pts(:,2)<=0 | pts(:,2)>imgsize(2),:)=[];

img2(sub2ind(imgsize,pts(:,1),pts(:,2)))=pts(:,3);