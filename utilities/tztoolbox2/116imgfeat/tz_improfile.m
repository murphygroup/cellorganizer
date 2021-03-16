function x = tz_improfile(img,s,t)
%TZ_IMPROFILE Get values along a line in an image
%   X = TZ_IMPROFILE(IMG,S,T) returns a vector of image intensities of
%   pixels on the line from position S to position T.
%   
%   See also TZ_SETIMGPITSPIXEL

%   10-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

pts = ml_getlinepts(s,t);
ind = sub2ind(size(img),pts(:,1),pts(:,2));
x = img(ind);