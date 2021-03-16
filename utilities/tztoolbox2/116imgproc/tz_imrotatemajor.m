function [img2,ra]=tz_imrotatemajor(img)
%TZ_IMROTATEMAJOR Rotate an image to line it up with major axis. Obsolete.
%   IMG2 = TZ_IMROTATEMAJOR(IMG) returns the image that is rotated from the
%   [image] IMG with its major axis.
%   
%   [IMG2,RA] = TZ_IMROTATEMAJOR(...) also returns the angle of the major
%   axis. RA is a 1x3 vector: [orientation ((-45,135]), flipud (0/1),
%   fliplr (0/1)].
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2=tz_imrotatemajor(img)
%
%OVERVIEW:
%   rotate an image to line it up with major axis
%PARAMETERS:
%   img - original image
%RETURN:
%   img2 - rotated image
%   ra - rotate angle, 1x3 
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

mom = tz_bwmoment(img);
% 
% theta = .5 * atan((mom.mu02 - mom.mu20)/2/mom.mu11)+sign(mom.mu11)*pi/4+pi/2;
% %theta=-theta;

ra=zeros(1,3);

a=imfeature(img,'Orientation');
ra(1)=a.Orientation;

img2=imrotate(img, -a.Orientation-90, 'bilinear', 'loose');

img3=img2;
img3(:,sum(img3,1)==0)=[];
img3(sum(img3,2)==0,:)=[];

cx=ml_imgmoments(img3,1,0)/ ml_imgmoments(img3,0,0);
cy=ml_imgmoments(img3,0,1)/ ml_imgmoments(img3,0,0);


if cy<size(img3,1)/2
    img2=flipud(img2);
    img2=fliplr(img2);
    ra(2)=1;
end

if cx<size(img3,2)/2
    img2=fliplr(img2);
    ra(3)=1;
end

