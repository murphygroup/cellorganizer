function pts2 = tz_homocoord(pts)
%TZ_HOMOCOORD Convert points to homogenuous coordinates.
%   PTS2 = TZ_HOMOCOORD(PTS) returns the normalized coordinates of the 
%   points PTS, which is a DxN matrix if there are N points in D 
%   dimensional space.
%
%   See also TZ_HOMOEXT

%Sep-09-2005 Initial write T. ZHAO

if(nargin < 1)
    error('Exactly one argument is required.');
end

%under construction
lastRow = pts(end,:);
infIndices = find(lastRow==0);
lastRow(infIndices) = 1;
pts2 = pts./lastRow(ones(size(pts,1),1),:);
pts2(end,:) = 1;
pts2(end,infIndices) = 0;