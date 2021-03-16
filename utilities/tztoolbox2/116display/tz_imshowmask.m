function tz_imshowmask(img,mask)
%TZ_IMSHOWMASK Show image with mask edge.
%   TZ_IMSHOWMASK(IMG,MASK) show the superimposed image of IMG and the
%   edge of MASK, which is also an image. IMG and MASK should have the 
%   same size.

%   ??-???-2004 Initial write T. Zhao
%   01-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

mask=bwperim(mask);

img(mask==1)=max(img(:));
imshow(img,[]);
