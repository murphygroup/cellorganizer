function tz_plotptser_2d(pts,ptstyle,lnstyle)
%TZ_PLOTPTSER_2D Plot 2D polyline.
%   TZ_PLOTPTSER_2D(PTS,POINTTSTYLE,LINESTYLE) plots a polyline composed of 
%   points PTS, which is an Nx2 matrix for N points. 
%   POINTSTYLE and LINESTYLE are parameters for TZ_PLOTGRAPH.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

V=pts;

%make edges
for i=1:size(pts,1)-1
    E(i,:)=[i i+1];
end

tz_plotgraph(V,E,[],ptstyle,lnstyle,[]);


