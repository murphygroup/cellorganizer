function shape2 = tz_mxp2mxs(shape)
%TZ_MXP2MXS Obsolete. See ML_MXP2MXS.
%   SHAPE2 = TZ_MXP2MXS(SHAPE) returns the medial axis representation of
%   the shape which is originally represented by spline shape SHAPE.
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_mxp2mxs','ml_mxp2mxs'));

if nargin < 1
    error('Exactly 1 argument is required')
end

x=(0:shape.length-1)/(shape.length-1);

shape2.format = 'mxs';
shape2.medaxis = [(1:shape.length)',spval(shape.spmedaxis,x)'];
shape2.width = spval(shape.spwidth,x);

