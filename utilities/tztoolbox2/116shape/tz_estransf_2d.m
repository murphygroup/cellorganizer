function A = tz_estransf_2d(pts1,pts2,order)
%TZ_ESTRANSF_2D Estimate transformation between two sets of points.
%   A = TZ_ESTRANSF_2D(PTS1,PTS2,ORDER) returns the estimated 
%   transformation matrix from PTS1 to PTS2. ORDER is the order of
%   transformation. Strictly speaking, A  can only be called 'a 
%   transformation matrix' when order is 1. See TZ_EXPANDPTS for more
%   details about ORDER.
%
%   See also TZ_AFFINEPTS

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University



npt=size(pts1,1);
qpts1=tz_expandpts(pts1,order);

A(:,1)=regress(pts2(:,1),qpts1);
A(:,2)=regress(pts2(:,2),qpts1);

if order>0
    A=[A;0 0];
end

A(:,3)=[zeros(size(A,1)-1,1); 1];


