function img = tz_mxp2img(shape)
%TZ_MXP2IMG Obsolete. See ML_MXP2IMG.
%   IMG = TZ_MXP2IMG(SHAPE) returns an image that contains the medial axis
%   spline shape SHAPE.
%   
%   See also TZ_IMG2MXP

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_mxp2img','ml_mxp2img'));

if nargin < 1
    error('Exactly 1 argument is required')
end

shape2 = tz_mxp2mxs(shape);
img = tz_mxs2img(shape2.medaxis,shape2.width);
