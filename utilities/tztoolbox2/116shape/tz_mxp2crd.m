function pts = tz_mxp2crd(shape)
%TZ_MXP2CRD Obsolete. See ML_MXP2CRD.
%   PTS = TZ_MXP2CRD(SHAPE) returns an array of points from the medial axis
%   spline shape SHAPE.
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU


error(tz_genmsg('of','tz_mxp2crd','ml_mxp2crd'));

if nargin < 1
    error('Exactly 1 argument is required')
end

shape2 = tz_mxp2mxs(shape);
pts = tz_mxs2crd(shape2.medaxis,shape2.width);
