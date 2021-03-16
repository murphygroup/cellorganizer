function pts=tz_getptset_2d(s,dists,angles)
%TZ_GETPTSET_2D Get a series of points along line segments.
%   PTS = TZ_GETPTSET_2D(S,DISTS,ANGLES) returns 2D points (Nx2 matrix)
%   along line segments with lengths DISTS and angles ANGLES if the
%   length of vector DISTS is N. ANGLES has the same size as DISTS. The
%   line segments are connect from front to end, and their lengths and
%   directions are determined by DISTS and ANGLES.

%   16-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

pts(1,:)=s;
cangle=tz_findangle_2d([1,0],[0,0],s);

curpt=s;

for i=1:length(dists)
    cangle=cangle+angles(i);
    pts(i+1,:)=tz_rotate_2d([dists(i),0],cangle);
    pts(i+1,:)=pts(i+1,:)+pts(i,:);
end