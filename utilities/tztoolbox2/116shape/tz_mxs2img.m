function img = tz_mxs2img(medaxis,width)
%TZ_MXS2IMG Obsolete. See ML_MXS2IMG.
%   IMG = TZ_MXS2IMG(MEDAXIS,WIDTH) returns an image that contains the
%   shape represented by medial axis MEDAXIS and width WIDTH.
%   
%   See also

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_mxs2img','ml_mxs2img'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

img = tz_crd2img(tz_mxs2crd(medaxis,width));
