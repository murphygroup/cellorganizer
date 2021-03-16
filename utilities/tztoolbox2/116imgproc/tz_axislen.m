function len = tz_axislen(pts,a)
%TZ_AXISLEN The length of object expansion at a certain angle.
%   LEN = TZ_AXISLEN(PTS,A) returns the length of axis of the object PTS
%   at a certain angle A, which has the unit radius. PTS is a 2-column
%   matrix for 2D points.

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(pts,2)~=2
    error('Please input 2D points')
end

theta=a*pi/180;
mvec=[cos(theta),sin(theta)];
mproj=pts*mvec';
len=max(mproj)-min(mproj);