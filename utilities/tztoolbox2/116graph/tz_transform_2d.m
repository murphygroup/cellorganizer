function pt=tz_transform_2d(orgpt,tm)
%TZ_TRANSFORM_2D Linear transformation of 2D points.
%   PT = TZ_TRANSFORM_2D(ORGPT,TM) returns the points transformed from
%   the 2-column matrix ORGPT and the 3x3 or 3x2 transformation matrix TM.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

pt=[orgpt,ones(size(orgpt,1),1)]*tm;

pt=pt(:,1:2);