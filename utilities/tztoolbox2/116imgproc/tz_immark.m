function img2 = tz_immark(img,pts)
%TZ_IMMARK Mark image with points.
%  TZ_IMMARK(IMG,PTS) shows the image IMG with marked points at positions
%  PTS. Notice that PTS is from the image coordinate system.

%   14-Sep-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

imshow(img,[]);
hold on
impts=tz_im2figcoord(pts);
plot(impts(1,:),impts(2,:),'rx');
hold off