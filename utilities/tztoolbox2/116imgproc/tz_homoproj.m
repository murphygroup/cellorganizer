function pts2 = tz_homoproj(pts,h)
%TZ_HOMOTRANSFM Homogeneous projection.
%   PTS2 = TZ_HOMOPROJ(PTS,H) returns the transformed points by applying
%   transformation matrix H on points PTS, which is a DxN matrix if there 
%   are N points in D dimensional space. H must have D+1 columns. PTS2 is
%   a (M-1)xN matrix if H has M rows.

%   11-Sep-2005 Initial write T. Zhao
%   12-Sep-2005 Modified T. Zhao
%       - change function name

if nargin < 2
    error('Exactly 2 arguments are required')
end

pts2 = tz_homocoord( h*tz_homoext(pts) );
pts2(end,:)=[];

