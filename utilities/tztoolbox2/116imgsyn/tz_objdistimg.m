function img = tz_objdistimg(pts,imgsize)
%TZ_OBJDISTIMG Distance map of edge for an object.
%   IMG = TZ_OBJDISTIMG(PTS,IMGSIZE) returns an image IMG with IMGSIZE.
%   IMG is the distance map of the edge of the object PTS. Any pixel of
%   the object will take the opposite sign.

%   12-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

objimg=tz_obj2img(pts,imgsize,{'2d','bn'});

edgeimg=bwperim(objimg,8);

img=bwdist(edgeimg);

img(objimg>0)=-img(objimg>0);
