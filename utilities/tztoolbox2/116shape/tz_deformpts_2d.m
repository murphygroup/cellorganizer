function pts2 = tz_deformpts_2d(pts,A)
%TZ_DEFORMPTS_2D Obsolete. See ML_DEFORMPTS_2D.
%   PTS2 = TZ_DEFORMPTS_2D(PTS,A) returns an array of 2d points which is
%   the tranformation of PTS according to transformation matrix A.
%   
%   See also

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_deformpts_2d','ml_deformpts_2d'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

npt=size(pts,1);

order=-(size(A,1)-1)/2;

qpts=tz_expandpts(pts,order);

pts2=qpts*A;
pts2=pts2(:,1:2);
