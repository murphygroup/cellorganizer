function shape = tz_obj2mxp(obj)
%TZ_OBJ2MXP Convert an object to a medial axis spline representation.
%   SHAPE = TZ_OBJ2MXP(OBJ) returns the structure of medial axis spline
%   shape from the [object] or [point array] OBJ. 
%
%   See also

%   13-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[medaxis,width] = tz_obj2mxs(obj);
shape = tz_mxs2mxp(medaxis,width);