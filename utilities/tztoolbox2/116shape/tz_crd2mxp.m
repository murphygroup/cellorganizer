function shape = tz_crd2mxp(pts)
%TZ_CRD2MXP Convert coordinate shape into medial axis spline shape.
%   SHAPE = TZ_CRD2MXP(PTS) returns a structure of 'mxp' shape from the
%   coordinates PTS.
%   
%   See also

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[medaxis,width] = tz_crd2mxs(pts);

shape = tz_mxs2mxp(medaxis,width);