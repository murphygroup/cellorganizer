function curve2 = tz_closecurve(curve)
%TZ_CLOSECURVE Obsolete. See ML_CLOSECURVE.
%   CURVE2 = TZ_CLOSECURVE(CURVE) returns a [closed curve] from the [curve]
%   CURVE.
%   
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_closecurve','ml_closecurve'));

if nargin < 1
    error('Exactly 1 argument is required')
end

curve2 = curve;

if any(curve(1,:)~=curve(end,:))
    curve2(end+1,:) = curve(1,:);
end