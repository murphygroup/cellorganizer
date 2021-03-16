function pts2 = tz_fig2imcoord(pts)
%TZ_FIG2IMCOORD Convert figure coordinates to image coordinates.
%   PTS = TZ_IM2FIGCOORD(PTS) returns a point array representing image
%   coordinates of the point array PTS, which is supposed to be figure
%   coordinates. PTS is an 2xN matrix, within which each column 
%   represents one point.
%
%   See also TZ_IM2FIGCOORDS

%Sep-09-2005 Initial write T. ZHAO

if(nargin < 1)
    error('Exactly one argument is required.');
end

if(size(pts,1) ~= 2)
    error('The input must have exactly two rows');
end

pts2 = flipud(pts);