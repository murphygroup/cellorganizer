function dists2 = tz_shiftdist(dists,angle)
%TZ_SHIFTDIST Obsolete. See ML_SHIFTDIST.
%   DISTS2 = TZ_SHIFTDIST(DISTS,ANGLE) returns a vector of hit point
%   distances that are shifted from DISTS by ANGLE. This function assumes
%   that DISTS is started from angle 0 and the step is one degree.
%   
%   See also

%   10-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_shiftdist','ml_shiftdist'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

dists2 = circshift(dists,[0 -round(angle)]);