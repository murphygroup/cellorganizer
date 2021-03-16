function [medaxis,width] = tz_crd2mxs(pts)
%TZ_CRD2MXS Convert coordinates into medial axis.
%   MEDAXIS = TZ_CRD2MXS(PTS) returns the medial axis of the shape
%   coordinates PTS.
%   
%   [MEDAXIS,WIDTH] = TZ_CRD2MXS(...) also returns the width.
%   
%   See also TZ_IMAXIS TZ_MXS2CRD

%   30-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

pts = tz_showpts_2d(pts,'ln',1);
objimg = tz_obj2img(pts,[]);

[imgaxis,medaxis,width] = tz_imaxis(objimg);