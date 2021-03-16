function shape = tz_img2mxp(img)
%TZ_IMG2MXP Convert an image to shape.
%   SHAPE = TZ_IMG2MXP(IMG) returns the structure of medial axis spline
%   shape from the image IMG. This function will take all pixels in IMG
%   above intensity 0 as objects.
%   
%   See also TZ_MXP2IMG TZ_IMG2MXS

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[medaxis,width] = tz_img2mxs(img);
shape = tz_mxs2mxp(medaxis,width);