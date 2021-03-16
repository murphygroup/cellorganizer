function pts2 = tz_transform2_2d(pts,A)
%TZ_TRANSFORM2_2D 2D transformation of points.
%   PTS2 = TZ_TRANSFORM2_2D(PTS,A) returns the transformed points by
%   applying the transformation matrix A on the points PTS. A is a 3rd
%   order transformation matrix.
%   
%   See also

%   07-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

npt=size(pts,1);
qpts=[pts,pts.^2,pts.^3,ones(npt,1)];
pts2=qpts*A;
pts2=pts2(:,1:2);