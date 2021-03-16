function pts2 = tz_homotransfm(pts,h)
%TZ_HOMOTRANSFM Homogeneous transformation.
%   PTS2 = TZ_HOMOTRANSFM(PTS,H) returns the transformed points by applying
%   transformation matrix H on points PTS, which is a DxN matrix if there 
%   are N points in D dimensional space. H must have D+1 columns. PTS2 is
%   a (M-1)xN matrix if H has M rows.

%   11-Sep-2005 Initial write T. Zhao

error(tz_genmsg('of','tz_homotransfm','tz_homoproj'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

pts2 = tz_homocoord( h*tz_homoext(pts) );
pts2(end,:)=[];

