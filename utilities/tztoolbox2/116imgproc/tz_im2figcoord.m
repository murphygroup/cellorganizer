function pts2 = tz_im2figcoord(pts)
%TZ_IM2FIGCOORD Convert image coordinates to figure coordinates.
%   PTS = TZ_IM2FIGCOORD(PTS) returns a point array representing figure
%   coordinates of the point array PTS, which is supposed to be image
%   coordinates. PTS is an 2xN matrix, within which each column 
%   represents one point.
%
%   See also TZ_FIG2IMCOORDS

%Sep-09-2005 Initial write T. ZHAO

if(nargin < 1)
    error('Exactly one argument is required.');
end

if(size(pts,1) ~= 2)
    error('The input must have exactly two rows');
end

pts2 = flipud(pts);