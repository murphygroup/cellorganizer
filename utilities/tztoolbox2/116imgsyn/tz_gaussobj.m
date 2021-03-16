function obj = tz_gaussobj(sigma)
%TZ_GAUSSOBJ Obsolete. See ML_GAUSSOBJ.
%   OBJ = TZ_GAUSSOBJ(SIGMA) returns an object that is extracted from a
%   2D Gaussian distribution which has covariance matrix SIGMA. The object
%   contains no less than 95% energy of the Gaussian distribution.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_gaussobj','ml_gaussobj'));

if nargin < 1
    error('Exactly 1 argument is required')
end

img = tz_gaussimg(sigma);

y = tz_wquantile(img(:),0.95);
img(img<y) = 0;

imageSize = size(img);

[r,c]=find(img>0);
obj = [r,c,img(sub2ind(imageSize,r,c))];