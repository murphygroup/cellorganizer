function [dists,d1,d2] = tz_diffhitpts(pts1,pts2)
%TZ_DIFFHITPTS Calculate difference between two sets point by point.
%   DISTS = TZ_DIFFHITPTS(PTS1,PTS2) returns the the difference between two
%   sets of points, PTS1 and PTS2. Both PTS1 and PTS2 should be obtained
%   from boundary points at the same set of angles. DISTS is an Nx1 vector
%   if PTS1 and PTS2 have size Nx2.
%   
%   [DISTS,D1,D2] = TZ_DIFFHITPTS(...) also returns the distance between
%   PTS1 and point [0,0] and the distance between PTS2 and the point [0,0].
%   
%   See also

%   11-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University



%find start points
tpts1=tz_normapts(pts1,[0 0]);
tpts2=tz_normapts(pts2,[0,0]);

d1=sqrt(sum(tpts1.^2,2));
d2=sqrt(sum(tpts2.^2,2));
dists=d1-d2;
