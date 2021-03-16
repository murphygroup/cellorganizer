function pts2 = tz_homoext(pts)
%TZ_HOMOEXT Extend points to homogenuous coordinates
%   PTS2 = TZ_HOMOEXT(PTS) adds a row of ones at the last row of PTS to
%   make them homogenuous coordinates. PTS is a DxN matrix if there are
%   N points in D dimensional space.
%
%   See also TZ_HOMOCOORD

%Sep-09-2005 Initial write T. ZHAO

if(nargin < 1)
    error('Exactly one argument is required.');
end

pts2 = [pts;ones(1,size(pts,2))];