function A=tz_affinepts(pts1,pts2)
%TZ_AFFINEPTS Estimate affine transformation between two sets of points.
%   A = TZ_AFFINEPTS(PTS1,PTS2) returns the affine matrix for
%   transformation from PTS1 to PTS2. The estimation is based on least
%   square error. PTS1 and PTS2 have two columns and they must have the
%   same size. The number of points must be greater than 2.
%   
%   See also

%   ??-???-????  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

npt=size(pts1,1);
qpts1=[pts1,ones(npt,1)];
A(:,1)=regress(pts2(:,1),qpts1);
A(:,2)=regress(pts2(:,2),qpts1);
A(:,3)=[0 0 1]';